# human fraction 3 only

# steps
search 1: in mascot daemon, against whole human proteome
search 2: in mascot daemon, against 2 proteins only
export csv with 'Ions score cut-off' = 10 and Significance threshold = 0.05

# get all the identified peptides (one list per search)

# mapping
map the peptides of two csv lists to the human proteome
export in html format, with the identified peptides red-highlighted
