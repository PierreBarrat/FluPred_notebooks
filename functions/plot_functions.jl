using Distributions

function plot_all_trajectories(ph; yearticks=false, data=false)
	p = plot(size=(900,600))
	for z in ph
		X,Y,tmp,α = Flu.frequency_series(z)
		if yearticks
			X = [year(x) + month(x) /12. for x in X]
		end
		for a in 1:size(Y,2)
			plot!(p, X, Y[:,a], label="$(α[a])", linewidth=2.5)
		end
	end
	return p
end

function trajectory_freqbin(traj, alphabins)
	freqtraj_cf = Dict()
	for (α,dα) in alphabins
	    freqtraj_cf[α] = Flu.frequency_condition(traj, α, dα=dα)
	end
	return freqtraj_cf
end

function bernoulli_estimator(x,n)
	naivemean = x/n
	estmean = (x+1/2) / (n+1)
	α = x + 1/2
	β = n - x + 1/2
	P = Distributions.Beta(α,β)
	# 
	lowerbound = 0.
	while cdf(P, lowerbound) < 0.05
		lowerbound += 0.01
	end
	# 
	higherbound = 1.
	while cdf(P, higherbound) > 0.95
		higherbound -= 0.01
	end
	# 
	return estmean, estmean-lowerbound, higherbound-estmean
end

pfix(traj) = sum([x.fixation==:fixed for x in traj]) / (sum([x.fixation==:fixed for x in traj]) + sum([x.fixation==:lost for x in traj]));
meanfreq(traj) = mean(t.freq[t.index[:active]] for t in traj)
function pfix_v_freq(ph, alphabins)
	trajectories = Flu.all_trajectories(ph, keep_unfinished=false)
	trajectories = Flu.previous_state_condition(trajectories, :lost)
	# Binning by frequency
	traj_fb = sort(trajectory_freqbin(trajectories, alphabins));
	# Keeping only trajectories that have a frequency backed by 50 strains at the time where it is binned. 
	for (k,v) in traj_fb
	    traj_fb[k] = Flu.population_size_condition(v, 20, mode=:active)
	end
	# 
	n = [length(traj_fb[x]) for x in keys(traj_fb)] # For error bars
	x = [count(t->t.fixation==:fixed, traj_fb[x]) for x in keys(traj_fb)]
	out = bernoulli_estimator.(x,n)
	# 
	xf = [meanfreq(traj_fb[x]) for x in keys(traj_fb)]
	# yf = [x[1] for x in out]; 
	yf = [count(t->t.fixation==:fixed, traj_fb[x])/length(traj_fb[x]) for x in keys(traj_fb)]
	errup = [x[3] for x in out]; errdown = [x[2] for x in out]
	# 
	err = (errup, errdown)
	return xf, yf, err
end
function pfix_v_freq_positivederivative(ph, alphabins)
	trajectories = Flu.all_trajectories(ph, keep_unfinished=false)
	trajectories = Flu.previous_state_condition(trajectories, :lost)
	# Binning by frequency
	traj_fb = sort(trajectory_freqbin(trajectories, alphabins));
	# Keeping only trajectories that have a frequency backed by 50 strains at the time where it is binned. 
	for (k,v) in traj_fb
	    traj_fb[k] = Flu.population_size_condition(v, 20, mode=:active)
	    traj_fb[k] = Flu.derivative_condition(traj_fb[k])
	end
	# 
	n = [length(traj_fb[x]) for x in keys(traj_fb)] # For error bars
	x = [count(t->t.fixation==:fixed, traj_fb[x]) for x in keys(traj_fb)]
	out = bernoulli_estimator.(x,n)
	# 
	xf = [meanfreq(traj_fb[x]) for x in keys(traj_fb)]
	yf = [x[1] for x in out]; 
	errup = [x[3] for x in out]; errdown = [x[2] for x in out]
	# 
	err = (errup, errdown)
	return xf, yf, err
end


# function plot_ph(ph)
# 	p = plot(size=(900,600))
# 	X,Y,tmp,ab = 