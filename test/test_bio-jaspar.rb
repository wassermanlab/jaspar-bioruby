require 'helper'
require 'bio-jaspar'
require 'bio'

class TestBioJaspar < Test::Unit::TestCase
  context 'JASPAR module' do
  	should "correctly read jaspar formatted file" do
  		f = File.open('test/data/jaspar-test.jaspar', "r") 
			motifs = Bio::Motifs.parse(f, "jaspar")
			f.close

			# Test first motif in the set
			corr_motifs_beg_counts = {
				"A" => [3.0, 21.0, 25.0, 0.0, 0.0, 24.0, 1.0, 0.0], 
				"C" => [13.0, 1.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0], 
				"T" => [5.0, 3.0, 0.0, 25.0, 20.0, 0.0, 24.0, 23.0], 
				"G" => [4.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 2.0]
			}
			assert_equal corr_motifs_beg_counts, motifs[0].counts
			assert_equal 8, motifs[0].length
			assert_equal "HAT5", motifs[0].name
			assert_equal "MA0008.1", motifs[0].matrix_id

			# Test the last motif in the set
			corr_motifs_end_counts = {
				"A" => [4.0, 5.0, 5.0, 3.0, 0.0, 0.0, 25.0, 26.0, 0.0, 0.0, 26.0, 0.0, 0.0, 17.0, 0.0, 5.0, 2.0, 0.0, 0.0], 
				"C" => [2.0, 3.0, 4.0, 8.0, 1.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 25.0, 1.0, 6.0, 3.0, 5.0], 
				"T" => [2.0, 3.0, 5.0, 5.0, 0.0, 20.0, 0.0, 0.0, 26.0, 0.0, 0.0, 26.0, 23.0, 0.0, 0.0, 7.0, 4.0, 3.0, 0.0], 
				"G" => [1.0, 3.0, 2.0, 2.0, 25.0, 0.0, 1.0, 0.0, 0.0, 26.0, 0.0, 0.0, 3.0, 9.0, 1.0, 4.0, 0.0, 4.0, 3.0]
			}
			assert_equal corr_motifs_end_counts, motifs[-1].counts
			assert_equal 19, motifs[-1].length
			assert_equal "ATHB9", motifs[-1].name
			assert_equal "MA0573.1", motifs[-1].matrix_id
  	end

  	should "correctly read pfm formatted file" do
  		f = File.open('test/data/jaspar-test.pfm', "r") 
			motif = Bio::Motifs.parse(f, "pfm")
			f.close

			corr_counts = {
				"A" => [3.0, 21.0, 25.0, 0.0, 0.0, 24.0, 1.0, 0.0], 
				"C" => [13.0, 1.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0], 
				"T" => [5.0, 3.0, 0.0, 25.0, 20.0, 0.0, 24.0, 23.0], 
				"G" => [4.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 2.0]
			}
			assert_equal 1, motif.length # should only read one motif
			assert_equal 8, motif[0].length
			assert_equal corr_counts, motif[0].counts
  	end

  	should "correctly read sites formatted file" do
  		f = File.open('test/data/jaspar-test.sites', "r") 
			motif = Bio::Motifs.parse(f, "sites")
			f.close

			corr_counts = {
				"A" => [15, 4, 41, 36, 7, 19, 3], 
				"C" => [11, 35, 1, 2, 29, 14, 22], 
				"T" => [7, 2, 0, 1, 1, 3, 3], 
				"G" => [10, 2, 1, 4, 6, 7, 15]
			}
			assert_equal 7, motif[0].length
			assert_equal corr_counts, motif[0].counts

			# Check the first and last sequence motifs
			assert_equal "ccaaccc", motif[0].instances[0].to_s
			assert_equal "gtatctc", motif[0].instances[-1].to_s
  	end

  end

  # Once the reads pass the test, load the test files and setup the test
	setup do
		f = File.open('test/data/jaspar-test.jaspar', "r") 
		@motifs = Bio::Motifs.parse(f, "jaspar")
		@motif = @motifs.first
		f.close
	end

	context 'JASPAR module' do
  	should "correctly convert Motifs into jaspar formatted string" do
  		corr_jaspar = ">MA0008.1 HAT5\nA [  3.00  21.00  25.00   0.00   0.00  24.00   1.00   0.00]\nC [ 13.00   1.00   0.00   0.00   5.00   0.00   0.00   0.00]\nG [  4.00   0.00   0.00   0.00   0.00   1.00   0.00   2.00]\nT [  5.00   3.00   0.00  25.00  20.00   0.00  24.00  23.00]\n>MA0027.1 En1\nA [  4.00   5.00   3.00   0.00   4.00   3.00   3.00   2.00   1.00   1.00   1.00]\nC [  1.00   2.00   0.00   0.00   0.00   0.00   0.00   1.00   3.00   4.00   6.00]\nG [  2.00   2.00   7.00   2.00   3.00   7.00   0.00   4.00   3.00   1.00   1.00]\nT [  3.00   1.00   0.00   8.00   3.00   0.00   7.00   3.00   3.00   4.00   2.00]\n>MA0046.1 HNF1A\nA [  5.00   1.00   1.00   1.00  20.00  16.00   1.00   8.00  14.00   2.00   0.00  13.00   8.00   5.00]\nC [  0.00   0.00   0.00   0.00   0.00   2.00   0.00   2.00   0.00   0.00   4.00   1.00   8.00  13.00]\nG [ 14.00  20.00   0.00   0.00   0.00   1.00   0.00   4.00   1.00   0.00   0.00   3.00   3.00   0.00]\nT [  2.00   0.00  20.00  20.00   1.00   2.00  20.00   7.00   6.00  19.00  17.00   4.00   2.00   3.00]\n"
  		jaspar = Bio::Jaspar.write(@motifs[0, 3], "jaspar")
  		assert_equal corr_jaspar, jaspar
  	end

  	should "correctly convert Motifs into pfm formatted string" do
  		corr_pfm = "  3.00  21.00  25.00   0.00   0.00  24.00   1.00   0.00\n 13.00   1.00   0.00   0.00   5.00   0.00   0.00   0.00\n  4.00   0.00   0.00   0.00   0.00   1.00   0.00   2.00\n  5.00   3.00   0.00  25.00  20.00   0.00  24.00  23.00\n"
  		pfm = Bio::Jaspar.write(@motifs, "pfm")
  		assert_equal corr_pfm, pfm
  	end

  	should "correctly calculate pseudocounts" do
  		corr_pc = {"A" => 1.25, "C" => 1.25, "T" => 1.25, "G" => 1.25}
  		pc = Bio::Jaspar.calculate_pseudocounts(@motifs[0])
  		assert_equal corr_pc, pc
  	end
  end

  context 'JASPAR Motif class' do
		should "return a correct length" do
			assert_equal 8, @motif.length
		end

	  should "return a correct consensus sequence" do
	  	assert_equal Bio::Sequence.auto("CAATTATT").to_s, @motif.consensus.to_s
	  end

	  should "return a correct anticonsensus sequence" do
	  	assert_equal Bio::Sequence.auto("AGGGGTGA").to_s, @motif.anticonsensus.to_s
	  end

	  should "return a correct degenerate consensus" do
	  	assert_equal Bio::Sequence.auto("CAATTATT").to_s, @motif.degenerate_consensus.to_s
	  end

	  should "return a correct reverse complement" do
	  	corr_rc_counts = {
	  		"A" => [23.0, 24.0, 0.0, 20.0, 25.0, 0.0, 3.0, 5.0], 
	  		"C" => [2.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 4.0], 
	  		"T" => [0.0, 1.0, 24.0, 0.0, 0.0, 25.0, 21.0, 3.0], 
	  		"G" => [0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 1.0, 13.0]
	  	}
	  	rc = @motif.reverse_complement
	  	assert_equal 0.13, rc.counts.gc_content
	  	assert_equal 8, rc.length
	  	assert_equal corr_rc_counts, rc.counts
		end

		should "return a correct mask" do
			assert_equal [1,1,1,1,1,1,1,1], @motif.mask
		end

		should "return correct pseudocounts" do
			corr_pc = {"A" => 0.0, "C" => 0.0, "T" => 0.0, "G" => 0.0}
			assert_equal corr_pc, @motif.pseudocounts
		end

		should "return a correct background" do
			corr_bg = {"A" => 0.25, "C" => 0.25, "T" => 0.25, "G" => 0.25}
			assert_equal corr_bg, @motif.background
		end

		should "return a correct pwm" do
			corr_pwm = {
				"A" => [0.12, 0.84, 1.0, 0.0, 0.0, 0.96, 0.04, 0.0],
				"C" => [0.52, 0.04, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0],
				"T" => [0.2, 0.12, 0.0, 1.0, 0.8, 0.0, 0.96, 0.92],
				"G" => [0.16, 0.0, 0.0, 0.0, 0.0, 0.04, 0.0, 0.08]
			}
			assert_equal corr_pwm, @motif.pwm
		end

		should "return a correct pssm" do
			corr_pssm = {
				"A" => [-1.0588936890535685, 1.7484612330040357, 2.0, -Float::INFINITY, -Float::INFINITY, 1.9411063109464317, -2.643856189774725, -Float::INFINITY], 
				"C" => [1.0565835283663676, -2.643856189774725, -Float::INFINITY, -Float::INFINITY, -0.3219280948873623, -Float::INFINITY, -Float::INFINITY, -Float::INFINITY], 
				"T" => [-0.3219280948873623, -1.0588936890535685, -Float::INFINITY, 2.0, 1.6780719051126378, -Float::INFINITY, 1.9411063109464317, 1.8797057662822885], 
				"G" => [-0.6438561897747247, -Float::INFINITY, -Float::INFINITY, -Float::INFINITY, -Float::INFINITY, -2.643856189774725, -Float::INFINITY, -1.6438561897747248]
			}
			assert_equal corr_pssm, @motif.pssm
		end

		should "correctly format motif in jaspar format" do
			corr_jaspar_str = ">MA0008.1 HAT5\nA [  3.00  21.00  25.00   0.00   0.00  24.00   1.00   0.00]\nC [ 13.00   1.00   0.00   0.00   5.00   0.00   0.00   0.00]\nG [  4.00   0.00   0.00   0.00   0.00   1.00   0.00   2.00]\nT [  5.00   3.00   0.00  25.00  20.00   0.00  24.00  23.00]\n"
			assert_equal corr_jaspar_str, @motif.format("jaspar")
		end

		should "correctly format motif in pfm format" do
			corr_pfm_str = "  3.00  21.00  25.00   0.00   0.00  24.00   1.00   0.00\n 13.00   1.00   0.00   0.00   5.00   0.00   0.00   0.00\n  4.00   0.00   0.00   0.00   0.00   1.00   0.00   2.00\n  5.00   3.00   0.00  25.00  20.00   0.00  24.00  23.00\n"
			assert_equal corr_pfm_str, @motif.format("pfm")
		end

  end

  context "matrix" do
  	setup do
			@motif2 = @motifs[1]
			@non_inf_dist = @motifs[15].pssm.distribution
		end

  	should "correctly return maximum possible score" do
  		assert_equal 14.245035054658192, @motif.pssm.max
  	end

  	should "correctly return the minimum possible score" do
  		assert_equal -Float::INFINITY, @motif.pssm.min
  	end

  	should "correctly refuse fraction gc content calculation on pssm" do
  		assert_raise do
  			@motif.pssm.gc_content
  		end
  	end

  	should "correctly calculate the mean" do
  		assert_equal 11.882147864914165, @motif.pssm.mean
  	end

  	should "correctly calculate the std" do
  		assert_equal 2.315187013634166, @motif.pssm.std
  	end

  	should "correctly calculates the PWM score for the given sequence" do
  		corr_res = [-Float::INFINITY, -Float::INFINITY, -Float::INFINITY, 4.7579989]
  		res = @motif.pssm.calculate(Bio::Sequence.auto("AGTTAATTAAG")).map{ |a| 
  			if a.infinite?
  				a
  			else
  				(a * (10 ** 7)).floor / (10.0 ** 7)
  			end 
  		}
  		assert_equal corr_res, res
  	end

  	should "correctly search and return the position of the hits with PWM higher than threshold" do
  		corr_hits = [[3, 4.7579989]]
  		hits = @motif.pssm.search(Bio::Sequence.auto("AGTTAATTAAG")).map{ |a, b| 
  			[a, (b * 10 ** 7).floor / (10.0 ** 7)]
  		}
  		assert_equal corr_hits, hits
  	end

  	should "correctly compare sequences using Pearson's correlation" do
  		corr_pearson = [0.024879199790793116, -10]
  		assert_equal corr_pearson, @motif.pssm.dist_pearson(@motif2.pssm)
  	end

  	should "correctly generate a distribution for non-infinite pssms" do
  		assert_equal -54.665224748002345, @non_inf_dist.min_score
  		assert_equal 15000, @non_inf_dist.n_points
  		assert_equal 77.64489601267111, @non_inf_dist.interval
  		assert_equal 0.005176671512278893, @non_inf_dist.step

  		corr_md_beg_100 = [3.2600762748340303e-26, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.4153180022070797e-26, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.61062211083769e-26, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.782556497068055e-26, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.075095343542538e-26, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.269147502758849e-26, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.410691430657807e-26, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.5132776385471104e-26, 0.0, 0.0, 0.0, 0.0, 4.620724355927226e-26, 0.0, 0.0, 0.0, 0.0, 0.0]
  		corr_md_end_100 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.006362547198123186, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.007703187540387484]
  		assert_equal corr_md_beg_100, @non_inf_dist.mo_density[0, 100]
  		assert_equal corr_md_end_100, @non_inf_dist.mo_density[@non_inf_dist.mo_density.length-100, 100]

  		corr_bd_beg_100 = [9.313225746154785e-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.313225746154785e-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.313225746154785e-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.313225746154785e-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.313225746154785e-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.313225746154785e-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.313225746154785e-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.313225746154785e-10, 0.0, 0.0, 0.0, 0.0, 9.313225746154785e-10, 0.0, 0.0, 0.0, 0.0, 0.0]
  		corr_bd_end_100 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.313225746154785e-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.313225746154785e-10]
  		assert_equal corr_bd_beg_100, @non_inf_dist.bg_density[0, 100]
  		assert_equal corr_bd_end_100, @non_inf_dist.bg_density[@non_inf_dist.bg_density.length-100, 100]
  	end

  	should "correctly calculate the threshold for false positive rate" do
  		assert_equal -11.00517721344216, @non_inf_dist.threshold_fpr(0.1)
  	end

  	should "correctly calculate the threshold for false negative rate" do
  		assert_equal 8.655821190193073, @non_inf_dist.threshold_fnr(0.1)
  	end

  	should "correctly calculate the balanced threshold" do
  		assert_equal 0.3058500408872149, @non_inf_dist.threshold_balanced()
  	end
  	
  	should "correctly calculate the patser threshold" do
  		assert_equal 11.435693792286841, @non_inf_dist.threshold_patser()
  	end

  end
end
