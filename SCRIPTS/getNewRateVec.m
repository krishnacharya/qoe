function newRateVec = getNewRateVec(userVec, weightVec, C_0, maxUsersVec)
	% C_0 is nominal channel rate (in bits per slot)
	dr_lhs = sum(userVec);
	newRateVec = zeros(1, length(maxUsersVec));
	newRateVec(userVec > 0) = weightVec(userVec > 0) * C_0; 
	newRateVec = newRateVec / dr_lhs;
end