clear
nmsg = 500000;
nrcv = 3;
m = rand(nmsg, 1);
p0 = 0.5;
p1 = 1 - p0;
m(m < p0) = 0;
m(m >= p0) = 1;
mu = [0 0 0];
%rho = [1 0 0; 0 1 0; 0 0 1];
rho = [1 0.4 0.3; 0.4 1 0.4; 0.3 0.4 1];
sigma = [1 0 0; 0 1 0; 0 0 1];
for i = 1:nrcv
    for j = 1:nrcv
        if i ~= j
            sigma(i, j) = rho(i, j) * sqrt(sigma(i, i)) * sqrt(sigma(j, j));
        end
    end
end
w = 2 * (rho(1,2)^2 + rho(1,3)^2 + rho(2,3)^2 - 2*rho(1,2)*rho(1,3)*rho(2,3) - 1);
a = [rho(2,3)^2-1, rho(1,3)^2-1, rho(1,2)^2-1, rho(1,2)-rho(1,3)*rho(2,3), rho(1,3)-rho(1,2)*rho(2,3), rho(2,3) - rho(1,2)*rho(1,3)];
N = mvnrnd(mu, sigma, nmsg);
E = 0.1:0.1:15;
s = zeros(nmsg, 1);
R = zeros(nmsg, nrcv);
m_hat = zeros(nmsg, 1);
pe = zeros(1, 150);
for i = 1:150
    s(m == 0) = sqrt(E(i));
    s(m == 1) = -sqrt(E(i));
    s = repmat(s(:,1),1,3);
    R = N + s;
    R = R';
    thresh = -0.25 * w * log(p0/p1) / sqrt(E(i));
    r1 = R(1,:);
    r2 = R(2,:);
    r3 = R(3,:);
    likely = r1*a(1) + r2*a(2) + r3*a(3) + (r1+r2)*a(4) + (r1+r3)*a(5) + (r2+r3)*a(6);
    %likely = zeros(nmsg, 1);
    %for j = 1:nmsg
    %    r = R(j,:);
    %    likely(j) = r(1)*a(1) + r(2)*a(2) + r(3)*a(3) + (r(1)+r(2))*a(4) + (r(1)+r(3))*a(5) + (r(2)+r(3))*a(6);
    %end
    m_hat(likely >= thresh) = 1;
    m_hat(likely < thresh) = 0;
    pe(i) = (length(find(m ~= m_hat))/nmsg);
end

plot(E,log(pe),'DisplayName','H-D');

