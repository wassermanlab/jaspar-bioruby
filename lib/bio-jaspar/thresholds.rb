#--
# = bio-jaspar/thresholds.rb
# 
# Copyright:: (C) 2015-2015 Jessica Lee
# License:: Ruby License
#++

module Bio # :nodoc:
	module Motifs
		# This is a direct import of Bio.motifs.threshold module in Biopython.
		# The following document contains excerpts from Bio.motifs.threshold in 
		# Biopython:
		# 
		# ---
		# 
		# A class representing approximate score distribution for a given motif
		# 
    # Utilizes a dynamic programming approch to calculate the distribution of
    # scores with a predefined precision. Provides a number of methods for calculating
    # thresholds for motif occurences.
		class ScoreDistribution
			attr_reader :min_score, :interval, :n_points, :ic, :step, :mo_density,
									:bg_density

			# Create a score distribution instance
			def initialize(motif = nil, precision = 10 ** 3, pssm = nil, background = nil)
				if pssm.nil?
					@min_score = [0.0, motif.min_score].min
					@interval = [0.0, motif.max_score].max - @min_score
					@n_points = precision * motif.length
					@ic = motif.ic
				else
					@min_score = [0.0, pssm.min].min
					@interval = [0.0, pssm.max].max - @min_score
					@n_points = precision * pssm.length
					@ic = pssm.mean(background)
				end

				@step = @interval / (@n_points - 1)
				@mo_density = [0.0] * @n_points
				@mo_density[-_index_diff(@min_score)] = 1.0
				@bg_density = [0.0] * @n_points
				@bg_density[-_index_diff(@min_score)] = 1.0

				if pssm.nil?
					motif.log_odds.zip(motif.pwm).each do |lo, mo|
						self.modify(lo, mo, motif.background)
					end
				else
					(0...pssm.length).each do |position|
						mo_new = [0.0] * @n_points
						bg_new = [0.0] * @n_points
						lo = pssm[(0...pssm.alphabet.letters.length).to_a, position]
						lo.each do |letter, score|
							bg = background[letter]
							mo = (2 ** pssm[letter, position]) * bg
							d = _index_diff(score)
							
							(0...@n_points).each do |i|
								mo_new[_add(i, d)] += @mo_density[i] * mo
								bg_new[_add(i, d)] += @bg_density[i] * bg
							end
						end
						@mo_density = mo_new
						@bg_density = bg_new
					end
				end
			end

			# Modify motif density and background density
			def modify(scores, mo_probs, bg_probs)
				mo_new = [0.0] * @n_points
				bg_new = [0.0] * @n_points
				scores.each do |k, v|
					d = _index_diff(v)
					(0...@n_points).each do |i|
						mo_new[_add(i, d)] += @mo_density[i] * mo_probs[k]
						bg_new[_add(i, d)] += @bg_density[i] * bg_probs[k]
					end
				end
				@mo_density = mo_new
				@bg_density = bg_new
			end

			# Approximate the log-odds threshold which makes the type I error (false positive rate).
			def threshold_fpr(fpr)
				i = @n_points
				prob = 0.0
				while prob < fpr
					i -= 1
					prob += @bg_density[i]
				end
				return @min_score + i * @step
			end

			# Approximate the log-odds threshold which makes the type II error (false negative rate).
			def threshold_fnr(fnr)
				i = -1
				prob = 0.0
				while prob < fnr
					i += 1
					prob += @mo_density[i]
				end
				return @min_score + i * @step
			end

			# Approximate the log-odds threshold which makes FNR equal to FPR times rate_proportion
			def threshold_balanced(rate_proportion = 1.0, return_rate = false)
				i = @n_points
				fpr = 0.0
				fnr = 1.0
				while fpr * rate_proportion < fnr
					i -= 1
					fpr += @bg_density[i]
					fnr -= @mo_density[i]
				end

				if return_rate
					return @min_score + i * @step, fpr
				else
					return @min_score + i * @step
				end
			end

			# Threshold selection mimicking the behaviour of patser (Hertz, Stormo 1999) software.
			#
      # It selects such a threshold that the log(fpr)=-ic(M)
      # 
      # Note: the actual patser software uses natural logarithms instead of log_2, so the numbers
      # are not directly comparable.
			def threshold_patser
				return self.threshold_fpr(2 ** -@ic)
			end

			private

			def _index_diff(x, y = 0.0)
				# Ruby doesn't have double slashes for integer division;
				# Also, python double slashes return the floor of the float
				# Hence, return the floor of the resultig float
				return ((x - y + 0.5 * @step) / @step).floor
			end

			def _add(i, j)
				return [0, [@n_points - 1, i + j].min].max
			end
		end
	end
end
