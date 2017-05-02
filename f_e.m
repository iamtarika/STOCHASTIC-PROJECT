clear
nmsg = 500000;
m = rand(1, nmsg);
p0 = 0.5;
p1 = 1 - p0;
m(m >= p0) = 0;
m(m < p0) = 1;
v1 = 1;
v2 = 2;
v3 = 4;
n1 = sqrt(v1)*randn(1, nmsg);
n2 = sqrt(v2)*randn(1, nmsg);
n3 = sqrt(v3)*randn(1, nmsg);
E = 0.1:0.1:15;
s = zeros(1, nmsg);
r1 = zeros(1, nmsg);
r2 = zeros(1, nmsg);
r3 = zeros(1, nmsg);
m_hat = zeros(1, nmsg);
pe = zeros(1, 150);
for i = 1:150
    s(m == 0) = sqrt(E(i));
    s(m == 1) = -sqrt(E(i));
    r1 = s + n1;
    r2 = s + n2;
    r3 = s + n3;
    thresh = -0.5*log(p0/p1)/sqrt(E(i));
    likely = r1/v1 + r2/v2 + r3/v3;
    m_hat(likely >= thresh) = 0;
    m_hat(likely < thresh) = 1;
    pe(i) = (length(find(m ~= m_hat))/nmsg);
end

y=zeros(1,150);
y=log(pe);
plot(E,y,'DisplayName','F-E');
