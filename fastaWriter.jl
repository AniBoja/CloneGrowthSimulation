using CSV

for arg in ARGS
    all_mutations = CSV.File("out/all_mutations_" * arg * ".txt", delim="\t")
    number_mutations = length(all_mutations) + 1
    cell_mutations = CSV.File("out/cell_mutations_" * arg * ".tsv", delim="\t")
    fasta_output = open("cell_mutations_" * arg * ".fasta", "a")
    for row in cell_mutations
        if (row[4])
            cell_id = row[1]
            write(fasta_output, string(">", cell_id, "\n"))
            parent_mutations = []
            new_mutations = []
            mutations_present = zeros(Int, number_mutations)
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
            seq = replace(replace(chop(string(mutations_present); head=1, tail=1), ',' => ""), ' ' => "")
            wrapped_seq = replace(seq, r"(.{60})" => s"\1\n")
            write(fasta_output, wrapped_seq*"\n")
        end
    end
    close(fasta_output)
end


