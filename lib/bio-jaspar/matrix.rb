#--
# = bio-jaspar/matrix.rb
# 
# Copyright:: (C) 2015-2015 Jessica Lee
# License:: Ruby License
# 
# Biological matrix classes
# 
# A direct import of Bio.motifs.matrix module in Biopython
#++

require 'bio'

module Bio # :nodoc:
	module Motifs
		# A class representing general position matrix
		# 
		# <i>A direct import of Bio.motifs.matrix module in Biopython</i>
		class GenericPositionMatrix < Hash
			attr_reader :length, :alphabet
			
			# Create a GenericPostionMatrix instance
			def initialize(alphabet, values)
				super()

				@alphabet = nil
				alphabet.letters.each do |letter|
					if @length.nil?
						@length = values[letter].length
					elsif @length != values[letter].length
						raise IndexError, "data has inconsistent lengths"
					end
					self[letter] = values[letter].to_a	
				end
				@alphabet = alphabet
				@_letters = @alphabet.letters.sort
			end

			# Return a string representation of the position matrix
			def to_s
				words = (0...@length).map { |i| "%6d" % i }
				line = "   " + words.join(" ")
				lines = [line]
				@_letters.each do |letter|
					words = self[letter].map { |value| "%6.2f" % value }
					line = ("%c: " % letter) + words.join(" ")
					lines << line
				end
				text = lines.join("\n") + "\n"
				return text
			end

			# Hash getter by key
			# 
			# Ruby has no tuples so a tuple is assumed to be an array.
			# Also, Ruby has no slices. All conditionals involving a slice in
			# biopython is ignored here. It is expected that users pre-slice the
			# indices using Range and step method:
			# (start...stop).step(stride)
			# and feed the resulting array into this getter method.
			def [](*key)
				if key.is_a? Array
					if key.length == 2
						key1, key2 = key
						# Slice conditional in biopython is ignored
						if key1.is_a? Integer
							letter1 = @_letters[key1]
							dim1 = 1
						# Tuple in python is considered as an array here
						elsif key1.is_a? Array
							letters1 = key1.map { |i| @_letters[i] }
							dim1 = 2
						elsif key1.is_a? String
							if key1.length == 1
								letter1 = key1
								dim1 = 1
							else
								raise KeyError, "#{key1}"
							end
						else
							raise KeyError, "Cannot understand key #{key1}"
						end

						# Again, slice conditional is ignored, but in the case of
						# key2, we assume that the indices are already pre-sliced and
						# given as an array
						if key2.is_a? Array
							indices2 = key2
							dim2 = 2
						elsif key2.is_a? Integer
							index2 = key2
							dim2 = 1
						else
							raise KeyError, "Cannot understand key #{key2}"
						end
						
						if dim1 == 1 && dim2 == 1
							return Hash.instance_method(:[]).bind(self).call(letter1)[index2]
						elsif dim1 == 1 && dim2 == 2
							values = Hash.instance_method(:[]).bind(self).call(letter1)
							# Ruby has no tuple; replace with Array
							return indices2.map { |idx2| values[idx2] }
						elsif dim1 == 2 && dim2 == 1
							d = {}
							letters1.each { |ltr1| 
								d[ltr1] = Hash.instance_method(:[]).bind(self).call(ltr1)[index2]
							}
							return d
						else
							d = {}
							letters1.each do |ltr1|
								values = Hash.instance_method(:[]).bind(self).call(ltr1)
								d[ltr1] = indices2.map { |idx2| values[idx2] }
							end
							if letters1.sort == @_letters
								return self.class.new(@alphabet, d)
							else
								return d
							end
						end

					elsif key.length == 1
						key = key[0]
					else
						raise KeyError, "keys should be 1- or 2-dimensional"
					end

				end
				
				# Slice is ignored; it is assumed that the key is pre-sliced
				if key.is_a? Integer
					letter = @_letters[key]
					dim = 1
				elsif key.is_a? Array
					letters = key.each{ |i| @_letters[i] }
					dim = 2
				elsif key.is_a? String
					if key.length == 1
						letter = key
						dim = 1
					else
						raise KeyError, key
					end
				else
					raise KeyError, "Cannot understand key #{key}"
				end

				if dim == 1
					return Hash.instance_method(:[]).bind(self).call(letter)
				elsif dim == 2
					d = {}
					letters.each { |ltr| 
						d[ltr] = Hash.instance_method(:[]).bind(self).call(ltr)
					}
					return d
				else
					raise RuntimeError, "Should not get here"
				end
									
			end

			# Return the consensus sequence (most probable pattern to be generated
			# from this motif)
			def consensus
				sequence = ""
				(0...@length).each do |i|
					maximum = -Float::INFINITY
					for letter in @alphabet.letters
						count = self[letter][i]
						
						if count > maximum
							maximum = count
							sequence_letter = letter
						end
					end
					
					sequence += sequence_letter
				end

				return Bio::Sequence.auto(sequence)
			end

			# Return the least probable pattern to be generated from this motif
			def anticonsensus
				sequence = ""
				(0...@length).each do |i|
					minimum = Float::INFINITY
					for letter in @alphabet.letters
						count = self[letter][i]
						
						if count < minimum
							minimum = count
							sequence_letter = letter
						end
					end
					
					sequence += sequence_letter
				end

				return Bio::Sequence.auto(sequence)
			end

      # Return degenerate consensus of a motif
      # 
      # Following the rules adapted from
      # 
      # \D. R. Cavener: "Comparison of the consensus sequence flanking
      # translational start sites in Drosophila and vertebrates."
      # Nucleic Acids Research 15(4): 1353-1361. (1987).
      # 
      # The same rules are used by TRANSFAC.
			def degenerate_consensus
        degenerate_nucleotide = {
            "A" => "A",
            "C" => "C",
            "G" => "G",
            "T" => "T",
            "AC" => "M",
            "AG" => "R",
            "AT" => "W",
            "CG" => "S",
            "CT" => "Y",
            "GT" => "K",
            "ACG" => "V",
            "ACT" => "H",
            "AGT" => "D",
            "CGT" => "B",
            "ACGT" => "N"
        }

        sequence = ""
        (0...@length).each do |i|
        	nucleotides = self.keys.sort { |a, b|
        		# Reverse sort
        		self[b][i] <=> self[a][i]
        	}
        	counts = nucleotides.map { |c| self[c][i] }
        	
        	# Follow the Cavener rules:
        	if (counts[0] >= counts[1..-1].inject(:+)) &&
        			(counts[0] >= 2 * counts[1])
        		key = nucleotides[0]
        	elsif (4 * counts[0, 2].inject(:+)) > (3 * counts.inject(:+))
        		key = nucleotides[0, 2].sort.join("")
        	elsif counts[3] == 0
        		key = nucleotides[0, 3].sort.join("")
        	else
        		key = "ACGT"
        	end

        	nucleotide = degenerate_nucleotide[key]
        	sequence += nucleotide
        end

        return Bio::Sequence.auto(sequence)
			end
			
			# Return the fraction GC content
			def gc_content
				alphabet = @alphabet
				gc_total = 0.0
				total = 0.0
				(0...@length).each do |i|
					alphabet.letters.each do |letter|
						if letter =~ /C|G/
							gc_total += self[letter][i]
						end
						total += self[letter][i]
					end
				end
				return gc_total / total
			end

			# Return the reverse complement of given position matrix
			def reverse_complement
				values = {}
				values["A"] = self["T"].reverse
				values["T"] = self["A"].reverse
				values["G"] = self["C"].reverse
				values["C"] = self["G"].reverse
				alphabet = @alphabet
				return self.class.new(alphabet, values)
			end

		end

		# A class representing a frequency position matrix
		# 
		# <i>A direct import of Bio.motifs.matrix module in Biopython</i>
		class FrequencyPositionMatrix < GenericPositionMatrix
			# Create and return a position-weight matrix by normalizing the counts matrix.
			#
			# If pseudocounts is None (default), no pseudocounts are added
			# to the counts.
			#
			# If pseudocounts is a number, it is added to the counts before
			# calculating the position-weight matrix.
			#
			# Alternatively, the pseudocounts can be a dictionary with a key
			# for each letter in the alphabet associated with the motif.
			def normalize(pseudocounts = nil)
				counts = {}
				if pseudocounts.nil?
					@alphabet.letters.each { |letter| 
						counts[letter] = [0.0] * @length
					}
				elsif pseudocounts.is_a? Hash
					@alphabet.letters.each { |letter| 
						counts[letter] = [pseudocounts[letter].to_f] * @length
					}
				else
					@alphabet.letters.each { |letter| 
						counts[letter] = [pseudocounts.to_f] * @length
					}
				end

				(0...@length).each{ |i|
					@alphabet.letters.each { |letter| 
						counts[letter][i] += self[letter][i]
					}
				}

				# Actual normalization is done in the PositionWeightMatrix initializer
				return PositionWeightMatrix.new(@alphabet, counts)
			end
		end

		# A class representing a position weight matrix
		# 
		# <i>A direct import of Bio.motifs.matrix module in Biopython</i>
		class PositionWeightMatrix < GenericPositionMatrix
			# Create an instance of position weight matrix
			def initialize(alphabet, counts)
				super(alphabet, counts)

				(0...@length).each do |i|
					total = alphabet.letters.map { |letter| self[letter][i].to_f }.inject(:+)
					alphabet.letters.each{ |letter| 
						self[letter][i] /= total
					}
				end

				alphabet.letters.each { |letter|
					self[letter] = self[letter].to_a
				}
			end

			# Returns the Position-Specific Scoring Matrix.
			#
			# The Position-Specific Scoring Matrix (PSSM) contains the log-odds
			# scores computed from the probability matrix and the background
			# probabilities. If the background is None, a uniform background
			# distribution is assumed.
			def log_odds(background = nil)
				values = {}
				alphabet = @alphabet

				if background.nil?
					background = Hash[@_letters.map { |l| [l, 1.0] }]
				else
					background = Hash[background]
				end

				total = background.values.inject(:+)
				alphabet.letters.each do |letter|
					background[letter] /= total
					values[letter] = []
				end

				(0...@length).each do |i|
					alphabet.letters.each do |letter|
						b = background[letter]
						if b > 0
							p = self[letter][i]
							if p > 0
								logodds = Math.log(p / b, 2)
							else
								logodds = -Float::INFINITY
							end
						else
							p = self[letter][i]
							if p > 0
								logodds = Float::INFINITY
							else
								logodds = Float::NAN
							end
						end

						values[letter] << logodds
					end
				end
				return PositionSpecificScoringMatrix.new(alphabet, values)
			end
		end

		# A class representing a position specific scoring matrix
		# 
		# <i>A direct import of Bio.motifs.matrix module in Biopython</i>
		class PositionSpecificScoringMatrix < GenericPositionMatrix
			# Returns the PWM score for a given sequence for all positions.
			#
			# Notes:
			#
			# - the sequence can only be a DNA sequence
			# - the search is performed only on one strand
			# - if the sequence and the motif have the same length, a single
			#   number is returned
			# - otherwise, the result is a one-dimensional list or numpy array (array
			#   for ruby)
			def calculate(sequence)
				if @alphabet.letters.join("") =~ /[^ACGT]/
					raise ArgumentError, "PSSM has wrong alphabet: #{@alphabet} - Use only with DNA motifs"
				end
				if sequence.illegal_bases.length > 0
					raise ArgumentError, "Sequence has wrong alphabet: #{sequence} - Use only with DNA sequences"
				end

				sequence = sequence.to_s
				m = @length
				n = sequence.length

				scores = []
				sequence = sequence.upcase
				(0...(n - m + 1)).each do |i|
					score = 0.0
					(0...m).each do |position|
						letter = sequence[i + position]
						begin
							score += self[letter][position]
						rescue KeyError
							score = Float::NAN
							break
						end
					end
					
					scores << score
				end

				if scores.length == 1
					return scores[0]
				else
					return scores
				end
			end

			# Find hits with PWM score above given threshold.
			#
      # A generator function, returning found hits in the given sequence
      # with the pwm score higher than the threshold.
			def search(sequence, threshold = 0.0, both = true)				
	      # sequence = sequence.to_s.upcase # No need to upcase
	      Enumerator.new do |enum|
		      n = sequence.length
		      m = @length

		      if both
		      	rc = reverse_complement
		      end
		      (0...(n - m + 1)).each do |position|
		      	# Bioruby sequences are 1-based (the end doesn't need +1 because
		      	# subseq() includes the end position)
		      	s = sequence.subseq((position + 1), (position + m)) 
		      	score = calculate(s)
		      	if score > threshold
		      		enum.yield([position, score])
		      	end

		      	if both
		      		score = rc.calculate(s)
			      	if score > threshold
		      			enum.yield([position - n, score])
		      		end
		      	end
		      end
		    end
			end

			# Maximal possible score for this motif
			#
	    # returns the score computed for the consensus sequence.
			def max
				score = 0.0
				letters = @_letters
				(0...@length).each do |position|
					score += letters.map { |letter| self[letter][position] }.max
				end
				return score
			end

			# Minimal possible score for this motif
			#
	    # returns the score computed for the anticonsensus sequence.
			def min
				score = 0.0
				letters = @_letters
				(0...@length).each { |position|
					score += letters.map { |letter| self[letter][position] }.min
				}
				return score
			end

			# Return fraction GC content
			def gc_content
				raise "Cannot compute the %GC composition of a PSSM"
			end

			# Return expected value of the score of a motif
			def mean(background = nil)
				if background.nil?
					background = Hash[@_letters.map { |l| [l, 1.0] }]
				else
					background = Hash[background]
				end
				total = background.values.inject(:+)

				@_letters.each { |letter|
					background[letter] /= total
				}

				sx = 0.0
				(0...@length).each do |i|
					@_letters.each do |letter|
						logodds = self[letter, i]
						if logodds.nan?
							next
						end
						if logodds.infinite? && logodds.infinite? < 0
							next
						end

						b = background[letter]
						p = b * (2 ** logodds)
						sx += p * logodds
					end
				end

				return sx
			end

			# Standard deviation of the score of a motif
			def std(background = nil)
				if background.nil?
					background = Hash[@_letters.map { |l| [l, 1.0] }]
				else
					background = Hash[background]
				end
				total = background.values.inject(:+)

				@_letters.each { |letter|
					background[letter] /= total
				}

				variance = 0.0
				(0...@length).each do |i|
					sx = 0.0
					sxx = 0.0
					@_letters.each do |letter|
						logodds = self[letter, i]
						if logodds.nan?
							next
						end
						if logodds.infinite? && logodds.infinite? < 0
							next
						end

						b = background[letter]
						p = b * (2 ** logodds)
						sx += p * logodds
						sxx += p * logodds * logodds					
					end
					sxx -= sx * sx
					variance += sxx
				end				
				variance = [variance, 0].max
				return Math.sqrt(variance)
			end

			# Return the similarity score based on pearson correlation for the 
			# given motif against self.
			#
	    # We use the Pearson's correlation of the respective probabilities.
			def dist_pearson(other)
				if @alphabet != other.alphabet
					raise ArgumentError, "Cannot compare motifs with different alphabets"
				end

				max_p = -2
				for offset in ((-@length + 1)...other.length).to_a
					if offset < 0
						p = self.dist_pearson_at(other, -offset)
					else
						p = other.dist_pearson_at(self, offset)
					end

					if max_p < p
						max_p = p
						max_o = -offset
					end
				end

				return (1 - max_p), max_o
			end

			# Return the similarity score based on pearson correlation for the given 
			# motif against self with given offset applied
			#
	    # We use the Pearson's correlation of the respective probabilities.
			def dist_pearson_at(other, offset)
				letters = @_letters
				sx = 0.0 # sum x
				sy = 0.0 # sum y
				sxx = 0.0 # sum x^2
				sxy = 0.0 # sum x dot y
				syy = 0.0
				norm = [@length, (offset + other.length)].max * letters.length

				(0...[(@length - offset), other.length].min).each do |pos|
					xi = letters.map { |letter| self[letter, (pos + offset)] }
					yi = letters.map { |letter| other[letter, pos] }

					sx += xi.inject(:+)
					sy += yi.inject(:+)
					sxx += xi.map { |x| x * x }.inject(:+)
					sxy += xi.zip(yi).map { |(x, y)| x * y }.inject(:+)
					syy += yi.map { |y| y * y }.inject(:+)
				end
				sx /= norm
				sy /= norm
				sxx /= norm
				sxy /= norm
				syy /= norm

				numerator = sxy - sx * sy
				denominator = Math.sqrt((sxx - sx * sx) * (syy - sy * sy))
				return numerator / denominator
			end

			# Calculate the distribution of the scores at the given precision
			def distribution(background = nil, precision = 10 ** 3)
				if background.nil?
					background = Hash[ @_letters.map { |l| [l, 1.0] } ]
				else
					background = Hash[background]
				end

				total = background.values.inject(:+)
				@_letters.each{ |letter| background[letter] /= total }

				return ScoreDistribution.new(nil, precision, self, background)
			end

		end

	end
end