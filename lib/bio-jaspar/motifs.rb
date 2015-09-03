#--
# = bio-jaspar/motifs.rb
# 
# Copyright:: (C) 2015-2015 Jessica Lee
# License:: Ruby License
#++

require 'bio'

module Bio # :nodoc:
	# == Tools for sequence motif analysis
	# 
	# All methods are directly imported from Bio.motifs module in Biopython.
	# The following document is based on the Bio.motifs document
	module Motifs

		# Create a motif instance based on sequence motifs
		def Motifs.create(instances, alphabet = nil)
			instances = Instances.new(instances, alphabet)
			return Motif.new(alphabet, instances, nil)
		end

		# Parses an output file of motif finding programs
		# 
		# === Arguments
		# +handle+:: file object
		# +format+:: currently supported input format
		# 
		# *NOTE:* only JASPAR formats are supported
		# 
    # Currently supported format:
    # 
    # * pfm:           JASPAR-style position-frequency matrix
    # * jaspar:        JASPAR-style multiple PFM format
    # * sites:         JASPAR-style sites file
		def Motifs.parse(handle, format)
			format = format.downcase
			if format =~ /^(pfm|sites|jaspar)$/
				record = Bio::Jaspar.read(handle, format)
				return record
			else
				raise ArgumentError, "Unknown format #{format}"
			end
		end

		# Read one motif from a handle using a specified file-format
		# 
		# === Arguments
		# +handle+:: file object
		# +format+:: currently supported input format
		# 
		# *NOTE:* only JASPAR formats are supported
		# 
		# For list of supported formats, go to Motifs.parse()
		def Motifs.read(handle, format)
			format = format.downcase
			motifs = parse(handle, format)
			if motifs.length == 0
				raise ArgumentError, "No motifs found in handle"
			end
			if motifs.length > 1
				raise ArgumentError, "More than one motif found in handle"
			end
			motif = motifs[0]
			return motif
		end

		# A class representing instances of sequence motifs
		# 
		# <i>A direct import of Bio.motifs module in Biopython</i>
		class Instances < Array
			attr_accessor :length
			attr_reader :alphabet 

			# Create an instance of sequence motif instances
			def initialize(instances = nil, alphabet = nil)
				super()

				if instances.nil?
					instances = []
				end
				@length = nil
				instances.each do |instance|
					if @length.nil?
						@length = instance.length
					elsif @length != instance.length
						message = "All instances should have the same length (#{instance.length} found, #{@length} expected)"
						raise ArgumentError, message	
					end

					begin
						a = instance.alphabet	
					rescue NoMethodError
						next
					end

					if alphabet.nil?
						alphabet = a
					elsif alphabet != a
						raise ArgumentError, "Alphabets are inconsistent"
					end
				end

				if alphabet.nil? || alphabet.letters.nil?
					alphabet = Alphabet.new.IUPAC_unambiguous_dna
				end

				instances.each do |instance|
					unless instance.is_a? Bio::Sequence
						sequence = instance.to_s
						instance = Bio::Sequence.auto(sequence)
					end

					self << instance
				end
				@alphabet = alphabet
			end

			# Return string representation of all sequence motifs
			def to_s
				text = ""
				self.each { |instance| 
					text += instance.to_s + "\n" 
				}
				return text
			end

			# Return the number of occurrences of nucleotides in each position 
			# of the motif instances (i.e. frequency position matrix)
			def count
				counts = {}
				@alphabet.letters.each { |letter|
					counts[letter] = [0.0] * @length
				}
				self.each do |instance|
					instance.to_s.upcase.each_char.each_with_index { |letter, position|
						counts[letter][position] += 1
					}
				end
				return counts
			end

			# A generator function - searches through the given sequence for
			# motif instances and returns found motif instances and their positions in
			# the given sequence
			# 
			# === Arguments
			# +sequence+:: Bio::Sequence object
			def search(sequence)
				Enumerator.new do |enum|
					(0...(sequence.length - @length + 1)).each do |pos|
						self.each do |instance|
							# Bioruby sequences are 1-based (the end doesn't need +1 because
		      		# subseq() includes the end position)
							if instance.to_s == sequence.subseq((pos + 1), (pos + @length)).to_s
								enum.yield([pos, instance])
								break
							end
						end
					end
				end
			end

			# Return reverse complements of all motif instances
			def reverse_complement
				instances = Instances.new(nil, @alphabet)
				instances.length = @length
				self.each do |instance|
					instance = instance.reverse_complement
					instances << instance
				end
				return instances
			end
		end

		# A class resembling Bio.Alphabet in Biopython
		# Only contains IUPAC unambiguous DNA
		class Alphabet
			attr_reader :letters, :size

			# Create an Alphabet instance
			def initialize
				@letters = nil
				@size = nil
			end

			# IUPAC unambiguous DNA nucleotides
			def IUPAC_unambiguous_dna
				@letters = ["G","A","T","C"]
				@size = 4

				return self
			end
		end

		# A sequence motif class
		# 
		# <i>A direct import of Bio.motifs module in Biopython</i>
		class Motif
			attr_accessor :length, :counts
			attr_reader :name, :alphabet, :instances

			# Create a motif instance
			def initialize(alphabet = nil, instances = nil, counts = nil)
				@name = ""
				if !counts.nil? && !instances.nil?
					raise ArgumentError, "Specify either instances or counts, don't specify both"
				elsif !counts.nil?
					if alphabet.nil?
						alphabet = Alphabet.new.IUPAC_unambiguous_dna
					end
					@instances = nil
					@counts = Bio::Motifs::FrequencyPositionMatrix.new(alphabet, counts)
					@length = @counts.length
				elsif !instances.nil?
					@instances = instances
					alphabet = @instances.alphabet
					counts = @instances.count
					@counts = Bio::Motifs::FrequencyPositionMatrix.new(alphabet, counts)
					@length = @counts.length
				else
					@counts = nil
					@instances = nil
					@length = nil
					if alphabet.nil?
						alphabet = Alphabet.new.IUPAC_unambiguous_dna
					end
				end
				@alphabet = alphabet

				# For the following attributes, we do not access the instance variables
				# directly; must call via setter methods so the default values are
				# properly set
				self.pseudocounts = nil
				self.background = nil
				self.mask = nil
			end

			# Getters and setters for mask, background, and pseudocounts begin

			# Get mask
			def mask
				return @_mask
			end

			# Set mask
			def mask=(mask)
				if @length.nil?
					@_mask = []
				elsif mask.nil?
					@_mask = [1] * @length
				elsif mask.length != @length
					raise ArgumentError, "The length (#{mask.length}) of the mask is inconsistent with the length (#{@length}) of the motif"
				elsif mask.is_a? String
					@_mask = []
					mask.each_char do |char|
						if char == "*"
							@_mask << 1
						elsif char == " "
							@_mask << 0
						else
							raise ArgumentError, "Mask should contain only '*' or ' ' and not a '#{char}'"
						end
					end
				else
					@_mask = mask.map do |c| 
						# Ruby doesn't have boolean type casting nor does it have 
						# bool to int conversion
						# Here we interpret integer 0 or nil as false
						# the rest as true
						if c.nil? || c == 0
							0
						else
							1
						end
					end
				end						
			end

			# Get pseudocounts
			def pseudocounts
				return @_pseudocounts
			end

			# Set pseudocounts
			def pseudocounts=(value)
				@_pseudocounts = {}
				if value.is_a? Hash
					@_pseudocounts = Hash[@alphabet.letters.map { |letter| [letter, value[letter]] }]
				else
					if value.nil?
						value = 0.0
					end
					@_pseudocounts = Hash[@alphabet.letters.map { |l| [l, value] }]
				end
			end

			# Get background
			def background
				return @_background
			end

			# Set background
			def background=(value)
				if value.is_a? Hash
					@_background = Hash[@alphabet.letters.map { |letter| [letter, value[letter]] }]
				elsif value.nil?
					@_background = Hash[@alphabet.letters.map { |l| [l, 1.0] }]
				else
					if @alphabet.letters.sort != ["A", "C", "G", "T"]
						raise "Setting the background to a single value only works for DNA motifs (in which case the value is interpreted as the GC content"
					end					
					@_background = {} if @_background.nil?
					@_background["A"] = (1.0 - value) / 2.0
					@_background["C"] = value / 2.0
					@_background["G"] = value / 2.0
					@_background["T"] = (1.0 - value) / 2.0
				end
				total = @_background.values().inject(:+)
				@alphabet.letters.each { |letter| 
					@_background[letter] /= total
				}
			end

			# Return position weight matrix of the given motif
			def pwm
				return @counts.normalize(@_pseudocounts)
			end

			# Return position specific scoring matrix of the given motif
			def pssm
				return self.pwm.log_odds(@_background)
			end

			# Return a string representation of all sequence motifs
			# 
			# === Arguments
			# +masked+:: if true, output sequences will be masked
			def to_s(masked = false)
				text = ""
				if !@instances.nil?
					text += @instances.to_s
				end

				if masked
					(0...@length).each do |i|
						if @_mask[i]
							text += "*"
						else
							text += " "
						end
					end
					text += "\n"
				end
				return text
			end

			# Return the length of the motif
			def length
				if @length.nil?
					return 0
				else
					return @length
				end
			end

			# Return the reverse complement of the motif. All the attributes of 
			# the current motif instance will be reverse complemented. e.g. 
			# frequency position matrix will be reverse complemented 
			def reverse_complement
				alphabet = @alphabet
				if !@instances.nil?
					instances = @instances.reverse_complement
					res = Motif.new(alphabet, instances)
				else
					res = Motif.new(alphabet, nil)
					counts = {}
					counts["A"] = self.counts["T"].reverse
					counts["T"] = self.counts["A"].reverse
					counts["G"] = self.counts["C"].reverse
					counts["C"] = self.counts["G"].reverse
					res.counts = FrequencyPositionMatrix.new(alphabet, counts)
					res.length = @length
				end				
				res.mask = @_mask.reverse if !@_mask.nil?
				return res
			end

			# Return the consensus sequence (most probable pattern to be generated
			# from this motif)
			def consensus
				return @counts.consensus
			end

			# Return the least probable pattern to be generated from this motif
			def anticonsensus
				return @counts.anticonsensus
			end

			# Generate degenerate consesnsus sequence
			# 
			# Following the rules adapted from:
      # \D. R. Cavener: "Comparison of the consensus sequence flanking
      # translational start sites in Drosophila and vertebrates."
      # Nucleic Acids Research 15(4): 1353-1361. (1987).
      # 
      # The same rules are used by TRANSFAC.
			def degenerate_consensus
				return @counts.degenerate_consensus
			end

			# TODO - skip since it uses transfac format
			# def weblogo(fname, format = "PNG", version = "2.8.2", **kwds)
			# end

			# Return a string output of the motif in the specified format
			# 
			# === Arguments
			# +format+:: output format (see currently supported format)
			# 
			# *NOTE:* Currently supports only JASPAR formats
			# 
			# Currently supported format:
			# 
			# * pfm: JASPAR single Position Frequency Matrix
			# * jaspar: JASPAR multiple Position Frequency Matrix
			# 
			def format(format)
				if format =~ /^(pfm|jaspar)$/
					motifs = [self]
					return Bio::Jaspar.write(motifs, format)
				else
					raise ArgumentError, "Unknown format type #{format}"
				end
			end
		end

		# Return a string representation of multiple motifs in a given format
		# 
		# === Arguments
		# +motifs+:: list of motif objects
		# +format+:: output format (see currently supported format)
		# 
		# *NOTE:* Currently supports only JASPAR formats
		# 
		# Currently supported format:
		# 
		# * pfm: JASPAR simple single Position Frequency Matrix
		# * jaspar: JASPAR multiple PFM format
		# 
		def Motifs.write(motifs, format)
			format = format.downcase
			if format =~ /^(pfm|jaspar)$/
				Bio::Jaspar.write(motifs, format)
			else
				raise ArgumentError, "Unknown format type #{format}"				
			end
		end

	end
end
