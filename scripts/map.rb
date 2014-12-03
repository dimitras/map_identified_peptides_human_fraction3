# ruby map.rb

require './csv_parser_for_protein_hits'
require './fasta_parser'

# io files
# csvp = CSVParserForProteinHits.open("../csvs/F003917_piped.csv", 100)
# fastap = FastaParser.open("../fasta/human_vasp_p2ry12.fasta")
# fileHtml = File.new("map_2.html", "w+")
#for all human proteome
csvp = CSVParserForProteinHits.open("../csvs/F003916_piped.csv", 100)
fastap = FastaParser.open("../fasta/HUMAN.fasta")
fileHtml = File.new("map_1.html", "w+")

fileHtml.puts "<HTML><BODY><font face=\"courier\" size=1>" 

csvp.each_protein do |protein|
	if fastap.entry_by_id(protein)
		prot_entry = fastap.entry_by_id(protein)
		prot_seq_highlight = []
		prot_seq_highlighted = ""
		csvp.protein_hits_per_protein(protein).each do |hit|
			for i in hit.pep_start.to_i-1..hit.pep_end.to_i-1
				prot_seq_highlight[i] = 1
			end
		end
		for i in 0..prot_entry.seq.length-1
			if prot_seq_highlight[i] == 1
				prot_seq_highlighted += "<font color=\"red\">#{prot_entry.seq[i]}</font>"
			else
				prot_seq_highlighted += prot_entry.seq[i]
			end
			if (i+1) % 60 == 0 #&& i != 0
				prot_seq_highlighted += "<br/>"
			end
		end
		fileHtml.puts prot_seq_highlighted + "<br/>"
	end
end

fileHtml.puts "</font></BODY></HTML>"
fileHtml.close()