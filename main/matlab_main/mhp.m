function [asdf_out, wtMat, prefFRs, dlys, ino_dlys] = mhp(N, p_init, asdf_in, p_inp, alph, bet, low_nu, sig)

	global dt = 0.1;
	adj_mat = rand(N) < p_init;

	adj_mat(1:(N+1):N^2) = 0;

	lininds = find(adj_mat(:));

	th = 1.0;

	mx_dly = 5;

	exc_sc = ones(N, 1);

	inpMat = rand(asdf_in{end-2}, N);
	inpMat = 0.1 .* inpMat .* (inpMat < p_init);
	inpinds = find(inpMat(:));
	inp_rast = asdfToRaster(asdf_in, 1, 'row');
	lastSpkTs = zeros(1,N);

	ref_e = 3;
	ref_i = 2;

	dlys = randi(5, length(lininds));
	wtMats = cell(mx_dly,1);

	for ii=1:mx_dly
		dly_inds = lininds(dlys == ii);

		if (length(dly_inds) / (N*(N-1))) < 0.2
			wtMats{ii} = sparse(N,N);
		else
			wtMats{ii} = zeros(N);
		end
		wtMats{ii}(dly_inds) = 0.01.*rand(1, length(dly_inds));
		wtMats{ii}(1:uint32(0.2*n), :) = wtMats{ii}(1:uint32(0.2*n), :) .* -1; % inhibitory
	end

	frs = zeros(1, N);
	actV = zeros(N, mx_dly);
	t_ptr = 1;
	recNet = zeros(1,N);

	ref_timer = zeros(N, 1);
	scale_fac_raw = zeros(1,N);
	scale_fac = ones(1,N);

	global frMat = adj_mat;
	global frValp = zeros(N,1);
	global frVald = zeros(N,1);
	global denom = log(1+alph);
	global inDegs = sum(adj_mat~=0);
	global [rInds, cInds] = find(adj_mat);

	function pfr_chs = meta(pfrs, frs, adj_mat, alph, bet, low_nu)

		frMat = frMat ./ frMat;
		frMat = frMat .* frs; % implicit singleton expansion

		frValp = exp(-pfrs./(low_nu * bet));
		frVald = log(1 + (alph * pfrs)/low_nu)/denom;

		for ii=1:N
			frMat(rInds(cInds==ii), ii) = (nonzeros(frs(ii) - frMat(rInds(cInds==ii), ii))  / pfrs(ii);
			frMat(frMat(:,ii)<0, ii) = -exp(frMat(frMat(:,ii)<0, ii)) * frVald(ii);
			frMat(frMat(:,ii)>0, ii) = exp(-frMat(frMat(:,ii)>0, ii)) * frValp(ii);
		end

		pfrs = pfrs + dt * (sum(frMat) ./ inDegs) .* (1+sig*randn(1,N));




	end

	ei = (1:N) > (N*0.2);

	for ii=1:asdf_in{end}(2)

		inNet = inp_rast(:, ii)' * inpMat;

		recNet = recNet * 0; % clear

		for kk=1:mx_dly
			ind = (t_ptr - kk) + 1;
			if ind < 1
				ind = ind + mx_dly;
			end
			netCh = actV(:, ind)' * wtMats{kk};
		end
		t_ptr = t_ptr + 1;
		if t_ptr > mx_dly
			t_ptr = 1;
		end

		active = ref_timer <= 0;
		ref_timer = ref_timer - 1;

		actV(:, t_ptr) = recNet > th & active;

		ref_timer(ei & actV(:, t_ptr)) = ref_e;
		ref_timer(~ei & actV(:, t_ptr)) = ref_i;

		lastSpkTs(actV(:, t_ptr)) = ii * dt;

		


	end

















end