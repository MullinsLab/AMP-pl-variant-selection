using NextGenSeqUtils, StatsBase

const IUPACbool = Dict{Char,Array{Bool,1}}(Dict())
IUPACbool['A']=[true,false,false,false]
IUPACbool['C']=[false,true,false,false]
IUPACbool['G']=[false,false,true,false]
IUPACbool['T']=[false,false,false,true]
IUPACbool['U']=[false,false,false,true]
IUPACbool['R']=[true,false,true,false]
IUPACbool['Y']=[false,true,false,true]
IUPACbool['S']=[false,true,true,false]
IUPACbool['W']=[true,false,false,true]
IUPACbool['K']=[false,false,true,true]
IUPACbool['M']=[true,true,false,false]
IUPACbool['B']=[false,true,true,true]
IUPACbool['D']=[true,false,true,true]
IUPACbool['H']=[true,true,false,true]
IUPACbool['V']=[true,true,true,false]
IUPACbool['N']=[true,true,true,true];

const nuc_dict = Dict('A' => 1,'C' => 2,'G' => 3,'T' => 4)

function resolve_base(c)
    uc = uppercase(c)
    if uc in keys(IUPACbool)
        return ['A','C','G','T'][rand((1:4)[IUPACbool[uc]])]
    else
        return uc
    end
end

function resolve_seq(s)
    return join(resolve_base.(collect(s)))
end

function string2nums(s)
    return [get(nuc_dict,c,0) for c in uppercase(s)]
end

#This is when your elements are stored as (size,AASEQ) tuples (size first, so its easy to sort).
#Can also be (size,AASEQ,AASEQ_name) tuples, with any number of things associated in later positions
function non_gap_tuple_hamming(s1,s2)
    return non_gap_hamming(s1[2],s2[2])
end

function next_frame_ct_gaps(s;ix=0)
    ct = 0
    frame =""
    for i in 1:length(s)
        ix+=1
        if s[i] != '-'
            ct+=1
            frame=frame*"$(s[i])"
            if ct == 3
            return frame,ix,s[i+1:end]
            end
        end
    end
    return "",ix,""
end

function find_first_stop_pos_with_gaps(s;ix=0)
    frame, ix, s_j = next_frame_ct_gaps(s;ix=ix)
    #println(frame)
    if frame in ["TAA","TAG","","TGA"]
        return ix
    else
        find_first_stop_pos_with_gaps(s_j;ix=ix)
    end
end

function extract_region(ali_seqs,ref_profile)
    sample_profile = seqs2profile(uppercase.(ali_seqs));
    ali_ref, ali_sample = profile_affine_align(ref_profile, sample_profile,profile_cost);
    matched_inds = [i[1][1] for i in ali_sample] .!= '!';
    region_inds = [i[1][1] for i in ali_ref] .!= '!';
    region_start,region_end = findfirst(region_inds[matched_inds]),findlast(region_inds[matched_inds]);
    extracted_seqs = [s[region_start:region_end] for s in ali_seqs]
    return extracted_seqs,(region_start,region_end)
end

function robust_translate(s)
    s = replace(s,"-"=>"")
    s = s[1:3*Int64(floor(length(s)/3))]
    translate_to_aa(s)
end

function find_first_stop(AAseq)
    s = collect(AAseq)
    if '*' in s
        return findfirst(s .== '*')
    else
        return length(AAseq)
    end
end

#Various distances.
function non_gap_hamming(s1,s2)
    l = min(length(s1),length(s2))
    diff = 0
    for i in 1:l
        if s1[i] != '-' && s2[i] != '-'
            if s1[i] != s2[i]
                diff += 1
            end
        end
    end
    return diff
end

function normalized_non_gap_hamming(s1,s2)
   l = min(length(s1),length(s2))
   diff = 0
   count = 0
   for i in 1:l
       if s1[i] != '-' && s2[i] != '-'
           count += 1
           if s1[i] != s2[i]
               diff += 1
           end
       end
   end
   return diff/count
end

function mafft_stop_preserving_align(AA_seqs)
    modified_seqs = [replace(s,"*"=>"B") for s in AA_seqs]
    AA_ali = mafft_align(modified_seqs)
    replaced_seqs = [replace(s,"B"=>"*") for s in AA_ali]
    return replaced_seqs
end

"""
    mafft_consensus{T<:BioSequence}(seqs::Vector{T}; kwargs...)
Julia wrapper for mafft.
"""
function mafft_align(seqs; kwargs...)
    mktempdir() do mydir
        seqfile = string(mydir, "/sequences.fasta")
        mafftout = string(mydir, "/mafft.fasta")
        write_fasta(string(mydir, "/sequences.fasta"), seqs)
        NextGenSeqUtils.mafft(seqfile, mafftout; kwargs...)
        aligned = Array{String}(read_fasta(mafftout))
        if  sum(uppercase.(degap.(seqs)) .!= uppercase.(degap.(aligned))) > 0
            @error "Aligned seqs don't match input order!"
        else
            return aligned
        end
    end
end

#Shorter version - doesn't do re-assignment
#Expects input data to be sorted
function pick_minor_by_dist_capped(dat,radius,dist_func; starting_set = [], max_return = 4)
    for rad in radius:50 #setting a cap here
        output_set = deepcopy(starting_set)
        if length(output_set) == 0
            push!(output_set,dat[1])
        end
        for d in dat
            if (minimum([dist_func(r,d) for r in output_set]) > rad) & (d[1] > 1)
                push!(output_set,d)
            end
        end
        println("Distance threshold: $(rad); Number of variants: $(length(output_set))...")
        if length(output_set) <= max_return
            return output_set
            break
        end
    end
end

function pick_RE_thresh_and_minority_variants(nt_collection,env_start,env_end,re_start,re_end;thresh=0.095,dist=8)
    AA_seqs = robust_translate.([s[env_start:env_end] for s in nt_collection])
    AA_seqs = getindex.(split.(AA_seqs,'*'),1);
    ali_AAs = mafft_stop_preserving_align(AA_seqs);
    variants = [(v,k) for (k,v) in collect(countmap(ali_AAs))]
    variants = sort(variants;rev=true)
    keeps = (getindex.(variants,1) ./ length(ali_AAs) .>= thresh) .& (getindex.(variants,1) .> 1)
    if length(variants[keeps]) == 0
        starting_set = []
    else
        starting_set = variants[keeps]
    end
    final_set = pick_minor_by_dist_capped(variants, dist, non_gap_tuple_hamming; starting_set = starting_set);
    synth_seqs, synth_names = [], []
    for (i,(size,variant)) in enumerate(final_set)
            #pull env by matching translation
            variant_members = [ali_AAs .== variant]
            variant_seqs = [s[re_start:re_end] for s in nt_collection[variant_members...]]
            most_freq = mode(variant_seqs)
            push!(synth_seqs, uppercase(most_freq))
            push!(synth_names, "RE_p00$(i)s Variant=00$(i) Frequency=$(round(size/length(AA_seqs); digits = 3)) Total#families=$(size)")
    end
    return synth_names, synth_seqs
end

env_seqs = read_fasta("panels/env_column_stripped_panel.fasta")
env_profile = seqs2profile(uppercase.(env_seqs))

re_seqs = read_fasta("panels/HIV1_COM_2017_5970-8795_DNA_stripped.fasta")
re_profile = seqs2profile(uppercase.(re_seqs))

f = snakemake.input[2]
f_intact = snakemake.input[1]*"/ofPossibleInterest/"*replace(snakemake.wildcards["sample"],"REN" => "env")*"_NT.fasta"
intact_names, intact_seqs = read_fasta_with_names(f_intact);

seqnames, ali_seqs = read_fasta_with_names(f);
ali_seqs = resolve_seq.(ali_seqs) #DISAMBIGUATE!

println("Processing alignment $(basename(f))...")

#extract env, RE sequences using profile alignment
env_seqs,(env_start,env_end) = extract_region(ali_seqs,env_profile);
re_seqs,(re_start,re_end) = extract_region(ali_seqs, re_profile);

println("RE start: $(re_start), RE end: $(re_end)...")
println("ENV start: $(env_start), ENV end: $(env_end)...")

keeps = [replace(s,"REN" => "env") in intact_names for s in seqnames]

#First function checker for alignment based start, but first occuring stop codon, no matter where it happens.
#from_start_of_env_AA = robust_translate.(a[env_start:end] for a in ali_seqs)
#function_checked = env_function_check_2.(from_start_of_env_AA)
#reason_names = seqnames .* " " .* [i[3] for i in function_checked]
#keeps = [i[1] for i in function_checked]
#translated = [i[2] for i in function_checked]
#AA_keeps_ali = mafft_stop_preserving_align(translated[keeps])
#write_fasta(snakemake.output[2],AA_keeps_ali,names = seqnames[keeps])
#write_fasta(snakemake.output[3],degap.(env_seqs[keeps]); names = seqnames[keeps])
#write_fasta(snakemake.output[4],degap.(env_seqs[.!keeps]); names = reason_names[.!keeps])

#pick AA frequency threshold and minority variants by distance
synth_names,synth_seqs=pick_RE_thresh_and_minority_variants(ali_seqs[keeps],env_start,length(ali_seqs[1]),re_start,length(ali_seqs[1]));
trimmed_synth_seqs = []
synth_names=snakemake.params["SID"].*"_".*synth_names
stop_pos = find_first_stop_pos_with_gaps.(s[env_start:end] for s in synth_seqs);
for (i,seq) in enumerate(synth_seqs)
    trim = seq[re_start:env_start + stop_pos[i] - 1]
    push!(trimmed_synth_seqs, trim)
end

write_fasta(snakemake.output[1],degap.(trimmed_synth_seqs);names=synth_names)
println("Done.")
