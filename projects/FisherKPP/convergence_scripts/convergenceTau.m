function convergenceTau

m = [10; 20; 40; 80; 160];

tau = [0.1; 0.05; 0.025; 0.0125; 0.00625];

eL2 = [0.0235131; 0.0120338; 0.00574872; 0.00250075; 0.00084534];

eL2 = sqrt(tau) .* eL2;

p1 = polyfit(log(m), log(eL2), 1);      % result: p1 = [-1.6862, -0.9034]
p2 = polyfit(log(tau), log(eL2), 1);    % result: p2 = [1.6862, -0.9034]

figure()
loglog(m, eL2, '-o')
hold on
loglog(m, exp(p1(1)*log(m) + p1(2)))
xlabel('Number of time steps')
ylabel('L2 error with respect to solution of smallest time step size');
legend('L2 error', 'slope: -1.7');

figure()
loglog(tau, eL2, '-o')
hold on
loglog(tau, exp(p2(1)*log(tau) + p2(2)))
xlabel('Time step size')
ylabel('L2 error with respect to solution of smallest time step size');
legend('L2 error', 'slope: 1.7');
