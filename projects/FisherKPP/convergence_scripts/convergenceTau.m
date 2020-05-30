function convergenceTau

m = [10; 20; 40; 80; 160];

tau = [0.1; 0.05; 0.025; 0.0125; 0.00625];

eL2 = [0.0235131; 0.0120338; 0.00574872; 0.00250075; 0.00084534];

eL2 = sqrt(tau) .* eL2;

p = polyfit(log(m), log(eL2), 1);       % result: p = [-1.6862, -0.9034]

figure()
loglog(m, eL2, '-o')
hold on
loglog(m, exp(p(1)*log(m) + p(2)))
xlabel('Number of time steps')
ylabel('L2 error with respect to solution of smallest time step size');
legend('L2 error', 'slope: -1.7');
