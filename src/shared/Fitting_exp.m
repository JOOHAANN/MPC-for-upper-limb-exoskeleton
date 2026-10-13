% 生成一些示例数据
x = (0:0.1:2)';
y = 2 * exp(0.5 * x) + randn(size(x)) * 0.1; % 加入一些噪声

% 定义指数模型
expModel = fittype('a*exp(b*x) + c*exp(d*x) + e', 'independent', 'x', 'coefficients', {'a', 'b', 'c', 'd', 'e'});

% 进行拟合
[fitresult, gof] = fit(x, y, expModel);

% extract the fitting results
a = fitresult.a;
b = fitresult.b;
c = fitresult.c;
d = fitresult.d;
e = fitresult.e;

% 绘制结果
figure;
plot(fitresult, x, y);
title('指数函数拟合');
xlabel('x');
ylabel('y');
legend('数据', '拟合曲线');
