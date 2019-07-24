function newGenes = recombAndMut(agentGenes, elites, survivors, toDie, ...
 tBRange, hidRange, recombP, mutationP, flipProb)
    
    N = size(agentGenes,1)-2;
    newGenes = zeros(size(agentGenes,1), toDie);
    mnHid = hidRange(1);
    mxHid = hidRange(2);
    mnTB = tBRange(1);
    mxTB = tBRange(2);
    
	for kk=1:toDie
		p = rand;
		if p < recombP % Time to make some babies

			% Randomly select 2 survivors
            coupling = survivors(randperm(length(survivors),2));

            randMix = rand(N+2,1)<0.5; % uniformly random gene slection
            newChromosome = zeros(N+2,1);

            % mix and match genes w/ parents
            newChromosome(randMix) = agentGenes(randMix, coupling(1));
            newChromosome(~randMix) = agentGenes(~randMix, coupling(2));
            
            % introduce a copying error
            mutation = rand(N+2,1)<flipProb;
            newChromosome(1:N) = xor((newChromosome(1:N)==1), mutation(1:N)); 
            if mutation(end-1)
                nH = newChromosome(end-1);
                newChromosome(end-1) = randi([max([mnHid, nH-15]), min([nH+15, mxHid])]);
            end
            if mutation(end)
            	tb = newChromosome(end);
                newChromosome(end) = randi([max([mnTB, tb-5]), min([mxTB, tb+5])]);
            end
            
            newGenes(:,kk) = newChromosome;

        elseif p >= recombP && p <= (recombP+mutationP) % just clone and mutate an elite
			    
			    toMut = elites(randi([1 length(elites)]));
                if rand < 0.5     % half the time mutate the input selection
                    mutation = rand(N+2,1) < flipProb*5;
                    newGenes(1:N, kk) = xor(agentGenes(1:N, toMut) == 1, mutation(1:N));
                else
                    newGenes(1:N,kk) = agentGenes(1:N, toMut);
                end
                
                % Always alter the number of hidden units and the amount of time looked back
                nH = agentGenes(end-1, toMut);
                newGenes(end-1, kk) = randi([max([mnHid, nH-15]), min([nH+15, mxHid])]);
                
                tB = agentGenes(end, toMut);
                newGenes(end, kk) = randi([max([mnTB, tB-5]), min([tB+5, mxTB])]);

		else % make an entirely new random agent
				newGenes(:, kk) = randomAgent(N+2, 0.02, .75, [mnHid mxHid], [mnTB mxTB]);
		end
	end
end