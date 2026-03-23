inputTbl = table( (datetime(2020,1,1):datetime(2025,1,1))',...
    randn(size((datetime(2020,1,1):datetime(2025,1,1))')),...
    randn(size((datetime(2020,1,1):datetime(2025,1,1))')), ...
    'VariableNames',{'tableTime','Asset1_Ret','Asset2_Ret'});
inputTbl = inputTbl(:,["Asset1_Ret","Asset2_Ret"]);
IF_LAGGED_IDV = true;
x0 = [0,0,0];
x_upper = [];
x_lower = [];
options = optimoptions("fmincon","Display","off");

rng(0,"twister")
theta = linspace(0.1,0.9,9);
tao = linspace(0.1,0.9,9)';

if(IF_LAGGED_IDV)
    [~,~,h] = ksdensity(inputTbl{1:end-1,1},"Bandwidth","plug-in");
else
    [~,~,h] = ksdensity(inputTbl{:,1},"Bandwidth","plug-in");
end

if(h<0.05)
    h = 0.05; % As recommended by Sim and Zhou (2015)
end

for i = 1:length(theta)
    for j = 1:length(tao)
        % disp([i,j])
        func_hdl = @(beta)qqr.loss(beta,inputTbl,theta(i),tao(j),h,IF_LAGGED_IDV);
        [solution,~] = fmincon(func_hdl,x0,[],[],[],[],x_lower,x_upper,...
            [],options);
        itcpt(i,j) = solution(1);
        slp_1(i,j) = solution(2);
        slp_2(i,j) = solution(3);
    end
end


% TODO: 统一刻度，排列堆叠, 图片以正方形为主
figure()
colormap("gray")
subplot(2,2,1)
contourf(theta,tao,itcpt',20, "EdgeColor","none","FaceColor","flat")
title('Intercept')
xlabel('quantile of dependent variable \theta')
ylabel('quantile of independent variable \tau')
colorbar()
subplot(2,2,2)
contourf(theta,tao,slp_1',20, "EdgeColor","none","FaceColor","flat")
title('Slope')
xlabel('quantile of dependent variable \theta')
ylabel('quantile of independent variable \tau')
colorbar()
