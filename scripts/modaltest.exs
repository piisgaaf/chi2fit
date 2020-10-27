alias Chi2fit.Fit, as: F
alias Chi2fit.Utilities, as: U

data1 = [{2,55},{4,40},{6,50},{8,25},{10,10},{12,5},{14,2}]
|> Enum.map(fn {x,count}->List.duplicate(x,count) end)
|> Enum.concat

data2 = [{2,30},{4,20},{6,50},{8,35},{10,27},{12,11},{14,20},{16,16},{18,2}]
|> Enum.map(fn {x,count}->List.duplicate(x,count) end)
|> Enum.concat

# The size of bins for grouping the data
binsize = 2

# Number of probes to use in the chi2 fit
_hdata1 = U.to_bins data1, {binsize,0}
hdata2 = U.to_bins data2, {binsize,0}

out = File.open!('single_xyz.dat', [:write])

# Single modal fit
IO.puts "========================================================================"
IO.puts "Single modal"
IO.puts "========================================================================"

#model = %Distribution.Weibull{}
#initial = [{1.655,1.67},{7.8,7.83}]
#options = [probes: 1_000_000, smoothing: false, model: :linear, saved?: true, surfacefile: out]
#result = {_,parameters2,_,saved} = F.chi2probe hdata2, initial, {Distribution.cdf(model), &F.nopenalties/2}, options
#U.display result

#options = [{:probes,saved}|options]
#result = {_,cov,parameters2,_} = F.chi2fit hdata2, {parameters2, Distribution.cdf(model), &F.nopenalties/2}, 1000, options
#U.display(hdata2,model,result,options)

## Bimodal fit
IO.puts "========================================================================"
IO.puts "Bimodal fit"
IO.puts "========================================================================"

#model = %Distribution.BiModal{distribs: [%Distribution.Weibull{},%Distribution.Weibull{}]}
#initial = [{0,1},{0.1,10},{0.1,10},{0.1,10},{0.1,10}]
#options = [probes: 100_000, smoothing: false, model: :linear, saved?: true]
#result = {_,parameters4,_,saved} = F.chi2probe hdata2, initial, {Distribution.cdf(model), &F.nopenalties/2}, options
#U.display result

#options = [{:probes,saved}|options]
#result = {_,cov,parameters4,_} = F.chi2fit hdata2, {parameters4, Distribution.cdf(model), &F.nopenalties/2}, 40, options
#U.display(hdata2,model,result,options)


## Trimodal fit
IO.puts "========================================================================"
IO.puts "Trimodal fit"
IO.puts "========================================================================"

penalty = fn _,[w1,w2|_] -> cond do
  w1 < 0 -> 1_000_000_000
  w1 > 1 -> 1_000_000_000
  w2 < 0 -> 1_000_000_000
  w2 > 1 -> 1_000_000_000
  true -> 0
  end
end

model = %Distribution.TriModal{distribs: List.duplicate(%Distribution.Weibull{},3)}
initial = [{0.29,0.32},{0.5,0.52},{1.15,1.3},{2.9,3.0},{0.79,0.8},{6.0,6.05},{4.9,5.0},{13.8,13.85}]
options = [probes: 1_000_000, smoothing: false, model: :linear, saved?: true, surfacefile: out]
result = {_,parameters5,_,saved} = F.chi2probe hdata2, initial, {Distribution.cdf(model), penalty}, options
U.display result

options = [probes: saved, smoothing: false, model: :linear, saved?: true, grid?: false]
result = {_,_cov,_parameters5,_} = F.chi2fit hdata2, {parameters5, Distribution.cdf(model), penalty}, 100, options
U.display hdata2,model,result,options
