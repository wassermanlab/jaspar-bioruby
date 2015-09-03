#--
# = bio-jaspar/jaspar.rb 
# 
# Copyright:: (C) 2015-2015 Jessica Lee
# License:: Ruby License
# 
# JASPAR 2014 module
# 
# A direct import of Bio.motifs.jaspar module in Biopython
#++

require 'bio'

module Bio # :nodoc:
	# == JASPAR 2014 module
	# 
	# Provides read access to a JASPAR5 formatted database.
	# 
	# This module is a direct import of Bio.motifs.jaspar module in Biopython.
	# The following document contains excerpts from Bio.motifs.jaspar module 
	# in Biopython.
	module Jaspar

		# Unambiguous DNA bases
		DNA = Bio::Motifs::Alphabet.new.IUPAC_unambiguous_dna

		# JASPAR OUTPUT specific DNA bases
		JASPAR_ORDERED_DNA_LETTERS = ["A","C","G","T"] # Jaspar requires specific order for printouts

		# A subclass of Bio::Motifs::Motif used to represent a JASPAR profile.
		# 
    # Additional metadata information are stored if available. The metadata
    # availability depends on the source of the JASPAR motif (a 'pfm' format
    # file, a 'jaspar' format file or a JASPAR database).
    # 
    # <i>A direct import of Bio.motifs.jaspar module in Biopython</i>
		class Motif < Bio::Motifs::Motif
			attr_accessor :matrix_id, :collection, :tf_class, :tf_family, :species,
										:tax_group, :acc, :data_type, :medline, :pazar_id, :comment

			# Construct a JASPAR Motif instance
			# 
			def initialize(matrix_id, name, opts = {})
				opts = {
					:alphabet => DNA, 
					:instances => nil, 
					:counts => nil, 
					:collection => nil, 
					:tf_class => nil, 
					:tf_family => nil, 
					:species => nil, 
					:tax_group => nil, 
					:acc => nil, 
					:data_type => nil, 
					:medline => nil, 
					:pazar_id => nil, 
					:comment => nil
				}.merge(opts)

				super(opts[:alphabet], opts[:instances], opts[:counts])

				@name = name
				@matrix_id = matrix_id
				@collection = opts[:collection]
				@tf_class = opts[:tf_class]
				@tf_family = opts[:tf_family]
				@species = opts[:species]
				@tax_group = opts[:tax_group]
				@acc = opts[:acc]
				@data_type = opts[:data_type]
				@medline = opts[:medline]
				@pazar_id = opts[:pazar_id]
				@comment = opts[:comment]
			end

			# Return the JASPAR base matrix ID
			def base_id
				base_id, _ = Jaspar.split_jaspar_id(@matrix_id)
				return base_id
			end

			# Return the JASPAR matrix version
			def version
				_, version = Jaspar.split_jaspar_id(@matrix_id)
				return version
			end

			# Return a string represention of the JASPAR profile.
			#
      # We choose to provide only the filled metadata information.
			def to_s
				tf_name_str = "TF name\t#{@name}\n"
				matrix_id_str = "Matrix ID\t#{@matrix_id}\n"
				the_string = tf_name_str + matrix_id_str

				if @collection
					collection_str = "Collection\t#{@collection}\n"
					the_string += collection_str
				end
        if @tf_class
          tf_class_str = "TF class\t#{@tf_class}\n"
          the_string += tf_class_str
        end
        if @tf_family
          tf_family_str = "TF family\t#{@tf_family}\n"
          the_string += tf_family_str
        end
        if @species
          species_str = "Species\t#{@species.join(",")}\n"
          the_string += species_str
        end
        if @tax_group
          tax_group_str = "Taxonomic group\t#{@tax_group}\n"
          the_string += tax_group_str
        end
        if @acc
          acc_str = "Accession\t#{@acc}\n"
          the_string += acc_str
        end
        if @data_type
          data_type_str = "Data type used\t#{@data_type}\n"
          the_string += data_type_str
        end
        if @medline
          medline_str = "Medline\t#{@medline}\n"
          the_string += medline_str
        end
        if @pazar_id
          pazar_id_str = "PAZAR ID\t#{@pazar_id}\n"
          the_string += pazar_id_str
        end
        if @comment
          comment_str = "Comments\t#{@comment}\n"
          the_string += comment_str
        end
        matrix_str = "Matrix:\n#{counts}\n\n"
        the_string += matrix_str
        return the_string
			end

			# Return the hash key corresponding to the JASPAR profile
			#
      # Note: We assume the unicity of matrix IDs
			def hash
				return @matrix_id.hash
			end

			# Compare two JASPAR motifs for equality. Two motifs are equal if their
			# matrix_ids match
			def ==(other)
				return @matrix_id == other.matrix_id
			end

		end

		# Represent a list of JASPAR motifs.
		#
		# <i>A direct import of Bio.motifs.jaspar module in Biopython</i>
    # 
    # ==== Attributes
    #
    # * +version+ - The JASPAR version used
		class Record < Array
			# Construct a record instance
			def initialize
				super()
				@version = nil
			end

			# Return a string of all JASPAR motifs in the list
			def to_s
				return self.map { |the_motif| the_motif.to_s }.join("\n")
			end

			# Return the list of matrices as a hash (ruby equivalent of dict)
			# of matrices
			def to_h
				dic = {}
				self.each { |motif|
					dic[motif.matrix_id] = motif
				}
				return dic
			end
		end

		# Read motif(s) from a file in one of several different JASPAR formats.
		# 
    # Return the record of PFM(s).
    # Call the appropriate routine based on the format passed
		def Jaspar.read(handle, format)
			format = format.downcase
			if format == "pfm"
				record = _read_pfm(handle)
				return record
			elsif format == "sites"
				record = _read_sites(handle)
				return record
			elsif format == "jaspar"
				record = _read_jaspar(handle)
				return record
			else
				raise ArgumentError, "Unknown JASPAR format #{format}"
			end
					
		end

		# Return the representation of motifs in "pfm" or "jaspar" format.
		def Jaspar.write(motifs, format)
			letters = JASPAR_ORDERED_DNA_LETTERS
			lines = []
			if format == "pfm"
				motif = motifs[0]
				counts = motif.counts
				letters.each do |letter|
					terms = counts[letter].map { |value| "%6.2f" % value }
					line = "#{terms.join(" ")}\n"
					lines << line
				end
			elsif format == "jaspar"
				motifs.each do |m|
					counts = m.counts
					line = ">#{m.matrix_id} #{m.name}\n"
					lines << line

					letters.each do |letter|
						terms = counts[letter].map { |value| "%6.2f" % value }
						line = "#{letter} [#{terms.join(" ")}]\n"
						lines << line
					end
				end
			else
				raise ArgumentError, "Unknown JASPAR format #{format}"
			end
				
			text = lines.join("")
			return text	
		end

		# Return pseudocounts of a given JASPAR motif
		def Jaspar.calculate_pseudocounts(motif)
			alphabet = motif.alphabet
			background = motif.background

			total = 0
			(0...motif.length).each do |i|
				total += alphabet.letters.map { |letter| motif.counts[letter][i].to_f }.inject(:+)
			end

			avg_nb_instances = total / motif.length
			sq_nb_instances = Math.sqrt(avg_nb_instances)

			if background
				background = Hash[background]
			else
				background = Hash[alphabet.letters.sort.map { |l| [l, 1.0] }]
			end

			total = background.values.inject(:+)
			pseudocounts = {}

			alphabet.letters.each do |letter|
				background[letter] /= total
				pseudocounts[letter] = sq_nb_instances * background[letter]
			end

			return pseudocounts
		end

		# Utility function to split a JASPAR matrix ID into its component.
		#
    # Components are base ID and version number, e.g. 'MA0047.2' is returned as
    # ('MA0047', 2).
		def Jaspar.split_jaspar_id(id)
			id_split = id.split(".")

			base_id = nil
			version = nil

			if id_split.length == 2
				base_id = id_split[0]
				version = id_split[1]
			else
				base_id = id
			end

			return base_id, version
		end

		# Private methods
		private

		# Read the motif from a JASPAR .pfm file (PRIVATE).
		def Jaspar._read_pfm(handle)
			alphabet = DNA
			counts = {}

			letters = JASPAR_ORDERED_DNA_LETTERS
			letters.zip(handle).each do |letter, line|
				words = line.split
				if words[0] == letter
					words = words[1..-1]
				end
				counts[letter] = words.map(&:to_f)
			end

			motif = Motif.new(nil, nil, :alphabet => alphabet, :counts => counts)
			motif.mask = "*" * motif.length
			record = Record.new
			record << motif

			return record
		end

		# Read the motif from JASPAR .sites file (PRIVATE).
		def Jaspar._read_sites(handle)
			alphabet = DNA
			instances = []

			handle_enum = handle.to_enum

			handle.each do |line|
				unless line.start_with?(">")
					break
				end

				line = handle_enum.next
				instance = ""
				line.strip.each_char do |c|
					if c == c.upcase
						instance += c
					end
				end
				instance = Bio::Sequence.auto(instance)
				instances << instance
			end

			instances = Bio::Motifs::Instances.new(instances, alphabet)
			motif = Motif.new(nil, nil, :alphabet => alphabet, :instances => instances)
			motif.mask = "*" * motif.length
			record = Record.new
			record << motif

			return record
		end

		# Read motifs from a JASPAR formatted file (PRIVATE).
		# 
    # Format is one or more records of the form, e.g.::
    # 
    #   - JASPAR 2010 matrix_only format::
    #
    #             >MA0001.1 AGL3
    #             A  [ 0  3 79 40 66 48 65 11 65  0 ]
    #             C  [94 75  4  3  1  2  5  2  3  3 ]
    #             G  [ 1  0  3  4  1  0  5  3 28 88 ]
    #             T  [ 2 19 11 50 29 47 22 81  1  6 ]
    #
    #   - JASPAR 2010-2014 PFMs format::
    #
    #             >MA0001.1 AGL3
    #             0	3	79	40	66	48	65	11	65	0
    #             94	75	4	3	1	2	5	2	3	3
    #             1	0	3	4	1	0	5	3	28	88
    #             2	19	11	50	29	47	22	81	1	6
    #             
		def Jaspar._read_jaspar(handle)
			alphabet = DNA
			counts = {}

			record = Record.new

			head_pat = /^>\s*(\S+)(\s+(\S+))?/
			row_pat_long = /\s*([ACGT])\s*\[\s*(.*)\s*\]/
			row_pat_short = /\s*(.+)\s*/

			identifier = nil
			name = nil
			row_count = 0
			nucleotides = ["A","C","G","T"]
			handle.each do |line|
				line = line.strip
				
				head_match = line.match(head_pat)
				row_match_long = line.match(row_pat_long)
				row_match_short = line.match(row_pat_short)

				if head_match
					identifier = head_match[1]
					if head_match[3]
						name = head_match[3]
					else
						name = identifier
					end
				elsif row_match_long
					letter, counts_str = row_match_long[1..2]
					words = counts_str.split
					counts[letter] = words.map(&:to_f)
					row_count += 1
					if row_count == 4
						record << Motif.new(identifier, 
						                    name, 
						                    :alphabet => alphabet, 
						                    :counts => counts)
						identifier = nil
						name = nil
						counts = {}
						row_count = 0
					end
				elsif row_match_short
					words = row_match_short[1].split
					counts[nucleotides[row_count]] = words.map(&:to_f)
					row_count += 1
					if row_count == 4
						record << Motif.new(identifier, 
						                    name, 
						                    :alphabet => alphabet, 
						                    :counts => counts)
						identifier = nil
						name = nil
						counts = {}
						row_count = 0
					end
				end		
			end

			return record
		end

		private_class_method :_read_pfm, :_read_sites, :_read_jaspar

	end
end
