using CSV

cell_mutations = CSV.File("out/cell_mutations_100.tsv", delim="\t")
fasta_output = open("cell_mutations.fasta", "a")
for row in cell_mutations
    cell_id = row[1]
    write(fasta_output, string(">Cell", cell_id, "\n"))
    parent_mutations = []
    new_mutations = []
    mutations_present = zeros(Int, 266)
    try
        parent_mutations = parse.(Int, split(chop(row[2]; head=1, tail=1), ','))
    catch ArgumentError
        parent_mutations = []
    end
    for mut in parent_mutations
        mutations_present[mut] = 1
    end
    try
        new_mutations = parse.(Int, split(chop(row[3]; head=1, tail=1), ',')) 
    catch ArgumentError
        new_mutations = []
    end
    for mut in new_mutations
        mutations_present[mut] = 1
    end
    write(fasta_output, replace(replace(chop(string(mutations_present); head=1, tail=1), ',' => ""), ' ' => "")*"\n")
end

close(fasta_output)