function ParameterRelSensitivity = ...
    LocalSensitivityAnalysis( DESystem, System_YFP_Output, Parameters, Y0,...
    Tol, odeSolver, options)
% LocalSensitivityAnalysus  Estimate Relative Local Sensitivity of DESystem
% to Parameters
%
%       DESystem        
%                       Function with signature DY = ( T, Y, P )
%                       where (T,Y) is a point in the solution and P is the 
%                       (unknown) parameter vector that needs to have
%                       sensitivity determined. (e.g. RNAi_DE_System or
%                       CRISPRi_DE_System
%       OutputSystem    
%                       Takes the output vector from a solver and returns
%                       the output value that is associated with the 
%                       measurement. (e.g. CRISPRi_System_YFP_Output)
%       Parameters
%                       Parameters to vary
%       Y0
%                       Initial Conditions for the DE (probably best if
%                       this aligns well with your measured data!)
%       Tol
%                       Percent by which we vary each parameter before
%                       measuring sensitivity (default: 0.01).
%       odeType
%                       If not given, uses ode45.
%   Returns:
%       ParameterRelSensitivity
%                       A vector associating a parameter (index) to 
%                       a given value representing sensitivity. Higher
%                       absolute values mean higher sensitivity.
%                       

    if nargin < 5
        Tol = 0.01;
    end
    if nargin < 6
        odeSolver = @ode45;
    end
    if nargin < 7
        options = odeset();
    end

    T = [0 ; 24*60*60*2];
    
    [~,Y1] = odeSolver( DESystem, T, Y0, options, Parameters );
    YFP1 = System_YFP_Output(Y1);
    
    ParameterRelSensitivity = zeros(size(Parameters));
    for p = 1:length(Parameters)
        Z = ones(size(Parameters));
        Z(p) = Tol +  1;
        
        DP = Parameters(p) * Tol;
        
        [~,Y2] = odeSolver( DESystem, T, Y0, options, Parameters .* Z );
        YFP2 = System_YFP_Output(Y2);
        
        ParameterRelSensitivity(p) = ...
            (( YFP2(end) - YFP1(end) ) / DP )...
            * ( Parameters(p) /  YFP1(end) );
    end
end
