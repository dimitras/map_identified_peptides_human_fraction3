class ProteinHit
	attr_accessor :prot_hit_num, :prot_acc, :prot_desc, :prot_score, :prot_mass, :prot_matches, :prot_matches_sig, :prot_sequences, :prot_sequences_sig, :prot_len, :pep_query, :pep_rank, :pep_isbold, :pep_isunique, :pep_exp_mz, :pep_exp_mr, :pep_exp_z, :pep_calc_mr, :pep_delta, :pep_start, :pep_end, :pep_miss, :pep_score, :pep_expect, :pep_res_before, :pep_seq, :pep_res_after, :pep_var_mod, :pep_var_mod_pos, :pep_num_match, :pep_scan_title

	def initialize(prot_hit_num, prot_acc, prot_desc, prot_score, prot_mass, prot_matches, prot_matches_sig, prot_sequences, prot_sequences_sig, prot_len, pep_query, pep_rank, pep_isbold, pep_isunique, pep_exp_mz, pep_exp_mr, pep_exp_z, pep_calc_mr, pep_delta, pep_start, pep_end, pep_miss, pep_score, pep_expect, pep_res_before, pep_seq, pep_res_after, pep_var_mod, pep_var_mod_pos, pep_num_match, pep_scan_title)
		@prot_hit_num = prot_hit_num
		@prot_acc = prot_acc
		@prot_desc = prot_desc.to_s
		@prot_score = prot_score
		@prot_mass = prot_mass
		@prot_matches = prot_matches
		@prot_matches_sig = prot_matches_sig
		@prot_sequences = prot_sequences
		@prot_sequences_sig = prot_sequences_sig
		@prot_len = prot_len
		@pep_query = pep_query
		@pep_rank = pep_rank
		@pep_isbold = pep_isbold
		@pep_isunique = pep_isunique
		@pep_exp_mz = pep_exp_mz
		@pep_exp_mr = pep_exp_mr
		@pep_exp_z = pep_exp_z
		@pep_calc_mr = pep_calc_mr
		@pep_delta = pep_delta
		@pep_start = pep_start
		@pep_end = pep_end
		@pep_miss = pep_miss
		@pep_score = pep_score.to_f
		@pep_expect = pep_expect.to_f
		@pep_res_before = pep_res_before
		@pep_seq = pep_seq.to_s
		@pep_res_after = pep_res_after
		@pep_var_mod = pep_var_mod
		@pep_var_mod_pos = pep_var_mod_pos
		@pep_num_match = pep_num_match
		@pep_scan_title = pep_scan_title
	end
	
	def to_csv()
		hit = '"' + [@prot_hit_num, @prot_acc, @prot_desc, @prot_score, @prot_mass, @prot_matches, @prot_matches_sig, @prot_sequences, @prot_sequences_sig, @prot_len, @pep_query, @pep_rank, @pep_isbold, @pep_isunique, @pep_exp_mz, @pep_exp_mr, @pep_exp_z, @pep_calc_mr, @pep_delta, @pep_start, @pep_end, @pep_miss, @pep_score, @pep_expect, @pep_res_before, @pep_seq, @pep_res_after, @pep_var_mod, @pep_var_mod_pos, @pep_num_match, @pep_scan_title].join('","') + '"'
		return hit
	end

end
	