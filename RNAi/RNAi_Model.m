function steadystate = RNAi_Model
% RNAi_Model Wrapper for RNAi modelling/graphing
%
% The wrapper calls many of the scripts within this directory, including:
%
%   RNAi_DE_System.m:            model function to be used with ode15s
%   RNAi_System_YFP_Output.m:    gets YFP concentration from a state vector        


	% Species
    titles = {'sRNA'; 'mRNA'; 'YFP'; 'Hfq-mRNA'; 'Hfq'; 'Hfq-sRNA Complex'; 'Hfq-mRNA-sRNA Complex'; 'mRNA-sRNA'};
	nEqns = 8;
    
	% Initial Conditions
    % ** Note: These were obtained from running a simulation to steady
    %          state with no production of sRNA within it. That's why they
    %          seem oddly specific.
	R0=zeros(nEqns,1);
    R0(1)= 0; % sRNA
    R0(2)= 2.0; % mRNA - same as CRISPR
    R0(3)= 20.0; % YFP
    R0(4)= 0; % Hfq-mRNA
    R0(5)= 0; % Hfq
    R0(6)= 0; % Hfq-sRNA Complex
	R0(7)= 0; % Hfq-sRNA-mRNA Complex
    R0(8)= 0; % bound sRNA-mRNA
    
    % Parameters & Explanation
    % * rate of productions - alpha
    % * rate of degradation - beta
	% * lowercase means RNA, uppercase mean protien product
    alpha_h = 1/600;            %* nM/s
    beta_h = 2.31*10^(-3);      % s^-1
    alpha_H = 0.9;              % s^-1
    beta_H = 6.42*10^(-5);      % s^-1 (degradation not provided?)
    beta_Hs = 6.42*10^(-5);     % s^-1 (degradation not provided?)
    beta_Hms = 6.42*10^(-5);    % s^-1 (degradation not provided?)
    alpha_s = 1/600;            %* nM/s
    beta_s = 2.31*10^(-3);      %* nM/s
    alpha_m = 1/600;            %* nM/s
    beta_m = 2.31*10^(-3);      %* s^-1
    alpha_M = 1/600;            %* s^-1
    beta_M = 6.42*10^(-5);      % s^-1
    k_1 = 10^(-3);              % (nM*s)^-1
    k_m1 = 0.7*10^(-4);         % s^-1 (the backward rate from Hfq-sRNA association was taken)
    %             (previously) to be zero. Hence this parameter didn't
    %             exist.
    k_2 = 3.5*10^(-3);          % (nM*s)^-1

    % degradosome is taken to be mass action
    % ** NOTE: This is where the major difference is.
    %          In a gram-negative cell there is a specialized 
    %          unit called a degradosome whose job it is to break
    %          down the Hfq-sRNA-mRNA complex. In gram-negative there isn't
    %          one, so the degradation is taken to be simply mass action.
    %          Without a degradosome, and with no data regarding the
    %          degradation of the Hfq-sRNA-mRNA complex, I'm not entirely
    %          sure what to do here.
    %    TL;DR: This should be mass action, and not Michaelis-Menton. This
    %           is where the confusion was.
    k_3 = 0.7*10^(-4);         % s^-1
    
    % New terms take into account binding of sRNA to mRNA w/o Hfq
    k_4 = k_2/10;               % (nM*s)^-1 (no value on wiki)
    beta_ms = 2.31*10^(-2);     % s^-1
   
    % Place all parameters in a vector
    parameters = [alpha_h beta_h alpha_H beta_H ...
                  alpha_s beta_s alpha_m beta_m alpha_M beta_M ...
                  k_1 k_m1 k_2 k_3 k_4 beta_ms beta_Hs beta_Hms];
    
	% ODE Solver and Plots
    RNAiODE=@RNAi_DE_System;
    options=odeset('Refine', 6);
	Tend = 100*60*60*1; % Time
    
    % Run simulation
    [t,R] = ode23s(RNAiODE, [0, Tend], R0, options, parameters);    
	figure(1);
    for i = 1:length(titles)
        subplot(2, 4, i);
        plot(t/3600, R(:,i), 'LineWidth', 2);
        title(titles(i));
        %axis tight
    end
    parameters(5) = 0; % Set alpha_S to zero = No sRNA activity
    [t,R] = ode23s(RNAiODE, [0, Tend], R0, options, parameters); 
    figure(2);
    for i = 1:length(titles)
        subplot(2, 4, i);
        plot(t/3600, R(:,i), 'LineWidth', 2);
        title(titles(i));
        axis tight
    end

    format long
    steadystate=R(end,:);

    % Local Sensitivity Analysis 
    RNAiRelSensitivity = ...
         LocalSensitivityAnalysis( RNAiODE, @RNAi_System_YFP_Output, ...
         parameters, R0, 0.01, 23 )
end