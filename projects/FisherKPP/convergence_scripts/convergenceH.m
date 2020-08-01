function convergenceH

N = [39; 125; 444; 1670];

meshsizes = [3.52544; 1.76272; 0.88136; 0.44068];

eL2 = [6.22193; 3.91591; 0.884941; 0.308175];

eL2 = sqrt(meshsizes) .* eL2;
 
p1 = polyfit(log(N), log(eL2), 1);          % result: p = [-1.1166, 6.7270]
p2 = polyfit(log(meshsizes), log(eL2), 1);  % result: p = [-2.0152, 0.1397]

figure()
loglog(N, eL2, '-o')
hold on
loglog(N, exp(p1(1)*log(N) + p1(2)))
xlabel('Number of degrees of freedom')
ylabel('L2 error with respect to solution on finest mesh');
legend('L2 error', 'slope: -1.12');

figure()
loglog(meshsizes, eL2, '-o')
hold on
loglog(meshsizes, exp(p2(1)*log(meshsizes) + p2(2)))
xlabel('Mesh Width h');
ylabel('L2 error with respect to solution on finest mesh');
legend('L2 error', 'slope: 2.02');
 