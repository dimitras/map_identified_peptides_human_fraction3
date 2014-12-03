require './protein_hit'

class CSVParserForProteinHits
	attr_accessor :filename, :cutoff

	def initialize(filename, cutoff)
		@filename = filename
		@cutoff = cutoff.to_f
		@filehandle = File.new(filename)
		@pep_index = {}
		create_pep_index
		@hits_index = Hash.new { |h,k| h[k] = [] }
		create_prot_index
	end

	def self.open(filename, cutoff)
		csvp = CSVParserForProteinHits.new(filename,cutoff)
		if block_given?
			csvp.each do |hit|
				yield hit
			end
		else
			return csvp
		end
	end


	def create_prot_index()
		in_protein_hits_table = false
		@filehandle.each do |line|
			line_pos = @filehandle.pos - line.length
			hit = line_parse(line)
			if hit.prot_hit_num.to_s == "prot_hit_num"
				in_protein_hits_table = true
				next
			elsif hit.prot_hit_num.to_s == "Peptide matches not assigned to protein hits"
				in_protein_hits_table = false
				break
			end
			if in_protein_hits_table == true
				@hits_index[hit.prot_acc] << line_pos
			end
		end
		@filehandle.rewind
	end

	def highest_from_cutoff_scored_hit_for_prot(protein)
		highest_scored_hit = nil
		hits = []
		hits = protein_hits_per_protein(protein)
		score = 0
		hits.each do |hit|
			if (hit.pep_expect <= @cutoff) && (hit.pep_score.to_f >= score)
				score = hit.pep_score
				if score >= @pep_score_cutoff
					highest_scored_hit = hit
				end
			end
		end
		return highest_scored_hit
	end

	def protein_hits_per_protein(protein)
		hits = []
		if !@hits_index.has_key?(protein)
			return hits
		end
		@hits_index[protein].each do |hit_pos|
			@filehandle.pos = hit_pos
			hit = line_parse(@filehandle.readline)
			hits << hit
		end
		return hits
	end

	def each_protein()
		@hits_index.each_key do |key|
			yield key
		end
	end


	def create_pep_index()
		in_protein_hits_table = false
		@filehandle.each do |line|
			line_pos = @filehandle.pos - line.length
			hit = line_parse(line)
			if hit.prot_hit_num.to_s == "prot_hit_num"
				in_protein_hits_table = true
				next
			end
			if in_protein_hits_table == true
				if !@pep_index.has_key?(hit.pep_seq)
					@pep_index[hit.pep_seq] = []
				end
				@pep_index[hit.pep_seq] << line_pos
			end
		end
		@filehandle.rewind
	end

	def highest_scored_hit_for_mod_pep(peptide, modification)
		highest_scored_hit = nil
		hits = []
		hits = protein_hits(peptide)
		score = 0
		hits.each do |hit|
			if (hit.pep_score.to_f >= score) && (hit.pep_var_mod.include? modification)
				score = hit.pep_score
				highest_scored_hit = hit	
			end
		end
		return highest_scored_hit
	end

	def highest_scored_hit_for_non_mod_pep(peptide, modification)
		highest_scored_hit = nil
		hits = []
		hits = protein_hits(peptide)
		score = 0
		hits.each do |hit|
			if (hit.pep_score.to_f >= score) && (!hit.pep_var_mod.include? modification)
				score = hit.pep_score
				highest_scored_hit = hit	
			end
		end
		return highest_scored_hit
	end

	def highest_scored_hit_for_pep(peptide)
		highest_scored_hit = nil
		hits = []
		hits = protein_hits(peptide)
		score = 0
		hits.each do |hit|
			if hit.pep_score.to_f >= score
				score = hit.pep_score
				highest_scored_hit = hit	
			end
		end
		return highest_scored_hit
	end

	def each()
		@pep_index.each_key do |key|
			hits = []
			@pep_index[key].each do |hit_pos|
				hit = hit_from_pos(hit_pos)
				hits << hit
			end
			yield hits
		end
	end
	
	def each_hit()
		@pep_index.each_key do |key|
			@pep_index[key].each do |hit_pos|
				yield hit_from_pos(hit_pos)
			end
		end
	end

	def each_peptide()
		@pep_index.each_key do |key|
			yield key
		end
	end

	def has_peptide(peptide)
		if @pep_index.has_key?(peptide)
			return true
		else
			return false
		end
	end

	def hit_from_pos(hit_pos)
		@filehandle.pos = hit_pos
		hit = line_parse(@filehandle.readline)
		# hit.rank = @rank_index[hit_pos.to_s]
		return hit
	end

	def protein_hits(peptide)
		hits = []
		if !@pep_index.has_key?(peptide)
			return hits
		end
		@pep_index[peptide].each do |hit_pos|
			hit = hit_from_pos(hit_pos)
			hits << hit
		end
		return hits
	end
	
	def line_parse(line)
		(prot_hit_num, prot_acc, prot_desc, prot_score, prot_mass, prot_matches, prot_matches_sig, prot_sequences, prot_sequences_sig, prot_len, pep_query, pep_rank, pep_isbold, pep_isunique, pep_exp_mz, pep_exp_mr, pep_exp_z, pep_calc_mr, pep_delta, pep_start, pep_end, pep_miss, pep_score, pep_expect, pep_res_before, pep_seq, pep_res_after, pep_var_mod, pep_var_mod_pos, pep_num_match, pep_scan_title) = line.chomp.split("|")
# 		prot_acc.gsub!("\"", "") if prot_acc
# 		prot_desc.gsub!("\"", "") if prot_desc
# 		pep_var_mod.gsub!("\"", "") if pep_var_mod
# 		pep_scan_title.gsub!("\"", "") if pep_scan_title
		return ProteinHit.new(prot_hit_num, prot_acc, prot_desc, prot_score, prot_mass, prot_matches, prot_matches_sig, prot_sequences, prot_sequences_sig, prot_len, pep_query, pep_rank, pep_isbold, pep_isunique, pep_exp_mz, pep_exp_mr, pep_exp_z, pep_calc_mr, pep_delta, pep_start, pep_end, pep_miss, pep_score, pep_expect, pep_res_before, pep_seq, pep_res_after, pep_var_mod, pep_var_mod_pos, pep_num_match, pep_scan_title)
	end

end
