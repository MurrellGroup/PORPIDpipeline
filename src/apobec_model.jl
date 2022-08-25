using Distributions
#All parameters in the log domain.
#All outputs in the log domain
const e = MathConstants.e

const nuc_dict = Dict('A' => 1,'C' => 2,'G' => 3,'T' => 4)
function string2nums(s)
    return [get(nuc_dict,c,0) for c in uppercase(s)]
end

function count_changes(s1,s2)
    num1 = string2nums(s1)
    num2 = string2nums(s2)
    counts = zeros(Int64,4,4)
    for i in 1:length(num1)
        if num1[i] > 0 && num2[i] > 0
            counts[num1[i],num2[i]] += 1
        end
    end
    return counts
end

function mu_mat(mu,ga; ti_tv = 4.5)
   #     a,c,g,t
   mat = [-(mu + ti_tv*mu + mu) mu ti_tv*mu mu;
          mu -(mu + mu + ti_tv*mu) mu ti_tv*mu;
          ti_tv*ga mu -(ti_tv*ga + mu + mu) mu;
          mu ti_tv*mu mu -(mu + ti_tv*mu + mu)]
   return mat
end

function Pmat(t,ga_mult)
    Q = mu_mat(1,ga_mult)
    P = exp(Q.*t)
    return P
end

function logLF(count_mat,t,ga_mult)
    P = Pmat(e^t,e^ga_mult)
    return sum(count_mat .* log.(P))
end

function logPrior(t,ga_mult)
    logpdf(Normal(-5,1.0),t) + log(0.99*pdf(Normal(0,0.1),ga_mult) + 0.01*pdf(Normal(0,1.0),ga_mult))
end

function APOBEC_debug(cons,ali_seq)
    count_mat = count_changes(cons,ali_seq)
    t_log_slices = collect(-12:0.5:-1);
    ga_mult_log_slices = collect(-1:0.05:5);

    posterior_grid = zeros(length(t_log_slices),length(ga_mult_log_slices))
    for i in 1:length(t_log_slices)
        for j in 1:length(ga_mult_log_slices)
            posterior_grid[i,j] = logLF(count_mat,t_log_slices[i],ga_mult_log_slices[j]) + logPrior(t_log_slices[i],ga_mult_log_slices[j])
        end
    end
    scaled = e.^(posterior_grid .- maximum(posterior_grid))
    scaled = scaled./sum(scaled);

    return count_mat,scaled, ga_mult_log_slices, sum(scaled,dims=1)[:]
end

function APOBEC(cons,ali_seq)
    count_mat = count_changes(cons,ali_seq)
    t_log_slices = collect(-12:0.5:-1);
    ga_mult_log_slices = collect(-1:0.05:5);

    posterior_grid = zeros(length(t_log_slices),length(ga_mult_log_slices))
    for i in 1:length(t_log_slices)
        for j in 1:length(ga_mult_log_slices)
            posterior_grid[i,j] = logLF(count_mat,t_log_slices[i],ga_mult_log_slices[j]) + logPrior(t_log_slices[i],ga_mult_log_slices[j])
        end
    end
    scaled = e.^(posterior_grid .- maximum(posterior_grid))
    scaled = scaled./sum(scaled);

    posterior_mean_ga_mult = e^sum(ga_mult_log_slices .* sum(scaled,dims=1)[:])
    posterior_ga_inflated = sum(sum(scaled,dims=1)[:][ga_mult_log_slices .> log(1.0)])
    posterior_mean_mutation_rate = e^sum(t_log_slices .* sum(scaled,dims=2)[:])

    return posterior_mean_ga_mult,posterior_ga_inflated,posterior_mean_mutation_rate,count_mat
end

function modelresults2table(model_results, seqnames)
    df = DataFrame([[el[i] for el in model_results] for i in 1:4],
                   [:posterior_mean_ga_mult,
                    :posterior_ga_inflated,
                    :posterior_mean_mutation_rate,
                    :count_mat]);
    df[!,:name] = seqnames;
    df[!,:ga_mutations] = [c[1] for c in get_mut_counts.(df[!,:count_mat])];
    df[!,:total_mutations] = [c[2] for c in get_mut_counts.(df[!,:count_mat])];
    df = df[!,[:name,:posterior_mean_ga_mult,:posterior_ga_inflated,
          :posterior_mean_mutation_rate,:ga_mutations,:total_mutations,:count_mat]];
    df = sort(df, :posterior_ga_inflated; rev = true)
    return df
end
