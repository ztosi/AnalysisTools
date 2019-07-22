function gene = randomAgent(N, mnProb, mxProb, hidRange, tbRange)
	gene = zeros(N,1);
	gene(1:N-2) = rand(N-2,1) < (rand*(mxProb-mnProb) + mnProb);
	gene(N-1) = randi(hidRange);
	gene(N) = randi(tbRange);
end