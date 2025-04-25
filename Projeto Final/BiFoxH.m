function out = BiFoxH(an1, alphan1, An1, ap1, alphap1, Ap1,...
    bq1, betaq1 ,Bq1, ...
    cn2, Cn2, cp2, Cp2, ...
    dm2, Dm2, dq2, Dq2, ...
    en3, En3, ep3, Ep3, ...
    fm3, Fm3, fq3, Fq3, ...
    x, y)


InftyValue = 50;
ZeroValue = 0.02;

%InftyValue = 5*max(abs([x, y]));
%InftyValue  = min([InftyValue, 100]);

%ZeroValue = min(abs([x, y]))/100;
% ZeroValue = max([ZeroValue, 0.002]);


Eval = 50000;
AbsTol = 1e-3;
RelTol = 1e-3;
Method = 'auto';
% Method = 'iterated';
%***** Integrand definition *****
F = @(s, t) ( GammaProd(1 - an1, alphan1, s, An1, t) .* GammaProd(dm2, -Dm2, s) .* GammaProd(1-cn2, Cn2, s) ...
    .* GammaProd(fm3, -Fm3, t) .* GammaProd(1 - en3, En3, t) .* (x.^s) .* (y.^t) ) ...
    ./ ( GammaProd(1 - bq1, betaq1, s, Bq1, t) .* GammaProd (ap1, -alphap1, s, -Ap1, t) ...
    .* GammaProd (1 - dq2, Dq2, s) .* GammaProd(cp2, -Cp2, s) .* GammaProd(1 - fq3, Fq3, t) ...
    .* GammaProd(ep3, -Ep3,t) );
%***** Contour definition *****

% Parameters:
p2 = length([Cn2 Cp2]);
q2= length([Dm2 Dq2]);
AstarS = sum(Cn2) - sum(Cp2) + sum(Dm2) - sum(Dq2);
DeltaS = sum([Dm2 Dq2]) - sum([Cn2 Cp2]);
deltaS = prod([Cn2 Cp2].^-[Cn2 Cp2]) * prod([Dm2 Dq2].^-[Dm2 Dq2]);
muS = sum([dm2 dq2]) - sum([cn2 cp2]) + (p2 - q2)/2;

% Conditions:
% contour Ls_(-infinity)
condition11S = (DeltaS>0)&& x~=0;
condition12S = (DeltaS==0) && abs(x)<deltaS && abs(x)>0;
condition13S = (DeltaS==0) && abs(x)==deltaS && real(muS) < -1;
condition1S = condition11S || condition12S || condition13S;
% contour Ls (+infinity)
condition21S = (DeltaS<0)&& x~=0;
condition22S = (DeltaS==0) && abs(x)>deltaS;
condition23S = (DeltaS==0) && abs (x) ==deltaS && real(muS) < -1;
condition2S = condition21S || condition22S || condition23S;
% Contour Ls_ (c+infinity),
condition3S  = AstarS >0 && abs(angle(x)) <AstarS*pi/2 && x~=0;
if(condition3S || (~condition1S && ~condition2S))
    %S contour Def
    epsilonS = ZeroValue;
    Sups = min(dm2./Dm2); Infs = max((cn2-1)./Cn2);
    if(isempty(Sups) && isempty(Infs))
        WpxS = 1;
    elseif(isempty(Sups) && ~isempty(Infs))
        WpxS = Infs + epsilonS;
    elseif(~isempty(Sups) && isempty(Infs))
        WpxS = Sups - epsilonS;
    else
        WpxS = (Sups + Infs)/2;% Sups<s<Infs
        %WpxS = (3*Sups + 2*Infs)/5;% Sups<s<Infs
    end
    
    infityS = InftyValue;
    % Conditions & Parapeters for t contour
    p3 = length([An1 Ap1 En3 Ep3]);
    q3 = length([Bq1 Fm3 Fq3]);
    AstarT = sum([An1 En3]) - sum([Ap1 Ep3]) + sum (Fm3) - sum([Bq1 Fq3]);
    DeltaT = sum([Bq1 Fm3 Fq3]) - sum( [An1 Ap1 En3 Ep3]);
    deltaT = prod([An1 Ap1 En3 Ep3].^-[An1 Ap1 En3 Ep3])*prod([Bq1 Fm3 Fq3].^-[Bq1 Fm3 Fq3]);
    muT=sum([(bq1 - WpxS * betaq1) fm3 fq3]) - sum([(an1-WpxS*alphan1) (ap1-WpxS*alphap1) en3 ep3 ])+(p3-q3)/2;
    
    % contour Lt_(-infinity)
    condition11T = (DeltaT>0) && y~=0;
    condition12T = (DeltaT==0) && abs(y) < deltaT && abs(y)>0;
    condition13T = (DeltaT==0) && abs(y) == deltaT && real(muT)<-1;
    condition1T = condition11T || condition12T || condition13T;
    % contour Lt_ (+infinity)
    condition21T = (DeltaT<0) && y~=0;
    condition22T = (DeltaT==0) && abs(y) > deltaT;
    condition23T = (DeltaT==0) && abs(y) == deltaT && real(muT)<-1;
    condition2T = condition21T || condition22T || condition23T;
    % Contour Lt (c+infinity),
    condition3T = AstarT > 0 && abs(angle(y)) < AstarT*pi/2 && y~=0;
    if (condition3T || (~condition1T && ~condition2T))
        epsilonT = ZeroValue;
        Supt = min(fm3./Fm3);
        Inft = max([((-1+an1-alphan1.*WpxS)./An1) ((en3-1)/En3)]);
        if(isempty(Supt) && isempty(Inft))
            WpxT=1;
        elseif(isempty (Supt) && ~isempty(Inft))
            WpxT = Inft +epsilonT;
        elseif (~isempty(Supt) && isempty(Inft))
            WpxT = Supt - epsilonT;
        else
            WpxT = (5*Supt + Inft)/6; %Supt<t<Inft
            
        end
        
        infityT = InftyValue;
        % quad2d settings
        Tol = 1e-6; MaxEval = Eval;
        %***** Bivariate Fox H *****
        X1 = WpxS-1i*infityS; X2 = WpxS+1i*infityS;
        Y1 = WpxT-1i*infityT; Y2 = WpxT+1i*infityT;
         
%         out = ((1/pi/2i)^2) * quad2d(F, X1, X2,Y1, Y2, 'AbsTol', Tol, 'RelTol', Tol, 'MaxFunEvals',...
%             MaxEval, 'Singular', true);
        out = ((1/pi/2i)^2) * integral2(F, X1, X2, Y1, Y2, 'AbsTol', AbsTol, 'RelTol', RelTol, 'Method', Method);
        return;
    end
    if (condition2T)
        epsilonT = 14*(y>10)+1;
        Supt = min(fm3./Fm3);
        Inft = max([((-1+an1-alphan1.*WpxS)./An1) ((en3-1)/En3)]);
        if(isempty (Supt) && isempty(Inft))
            WPxT=-1;
        elseif(isempty(Supt) && ~isempty(Inft))
            WPxT = Inft + 10*epsilonT;
        elseif(~isempty (Supt) && isempty (Inft))
            WPxT = Supt - epsilonT;
        else
            WPxT = (Supt + Inft)/2; % Supt<t<Infs
        end
        
        WayPointsT = [WPxT-1i*epsilonT WPxT+1i*epsilonT];
        infityT = -InftyValue;
        minT = min([((-1+an1-alphan1.*WpxS)./An1) ((en3-1)/En3)]);
        if(~isempty(minT))
            infityT = infityT + minT;
        end
    end
    
    if(condition1T)
        epsilonT = ZeroValue;
        Supt = min(fm3./Fm3); Inft = max([((-1+an1-alphan1.*WpxS)./An1) ((en3-1)/En3)]);
        if(isempty(Supt) && isempty(Inft))
            WPxT=1;
        elseif(isempty(Supt) && ~isempty(Inft))
            WPxT = Inft +10*epsilonT;
        elseif(~isempty (Supt) && isempty(Inft))
            WPxT = Supt -epsilonT;
        else
            WPxT = (Supt + Inft)/2;% Supt<t<Infs
        end
        
        WayPointsT = [WPxT-1i*epsilonT WPxT+1i*epsilonT];
        infityT = InftyValue;
        maxT= max(fm3./Fm3);
        if(~isempty(maxT))
            infityT = infityT + maxT;
        end
    end
    
    Tol = 1e-6; MaxEval = Eval;
    %***** Bivariate Fox H *****
    
    X1 = WpxS - 1i*infityS; X2 = WpxS + 1i*infityS;
    Y1 = infityT; Y2 = WayPointsT(1);
    
%     out1 = ((1/pi/2i)^2) * quad2d(F,X1, X2, Y1, Y2, 'AbsTol', Tol, 'RelTol', Tol, 'MaxFunEvals', ...
%         MaxEval, 'Singular', true);
    out1 = ((1/pi/2i)^2) * integral2(F, X1, X2, Y1, Y2, 'AbsTol', AbsTol, 'RelTol', RelTol, 'Method', Method);
    
    X1 = WpxS - 1i*infityS; X2 = WpxS + 1i*infityS;
    Y1 = WayPointsT(1); Y2 = WayPointsT(2);
%     out2 = ((1/pi/2i)^2) * quad2d(F, X1, X2, Y1, Y2, 'AbsTol', Tol, 'RelTol', Tol, 'MaxFunEvals', ...
%         MaxEval, 'Singular', true);
    out2 = ((1/pi/2i)^2) * integral2(F, X1, X2, Y1, Y2, 'AbsTol', AbsTol, 'RelTol', RelTol, 'Method', Method);
    
    X1 = WpxS - 1i*infityS; X2 = WpxS + 1i*infityS;
    Y1 = WayPointsT(2); Y2 = infityT;
%     out3 = ((1/pi/2i)^2) * quad2d(F, X1, X2, Y1, Y2, 'AbsTol', Tol, 'RelTol', Tol, 'MaxFunEvals', ...
%         MaxEval, 'Singular', true);
    out3 = ((1/pi/2i)^2) * integral2(F, X1, X2, Y1, Y2, 'AbsTol', AbsTol, 'RelTol', RelTol, 'Method', Method);
    
    out = out1 + out2 + out3;
    return;
end

%****** GammaProd Subfunction *******
    function output = GammaProd(p, x, X, y, Y)
        if(nargin==3)
            [pp, XX] = meshgrid(p, X);
            xx = meshgrid(x, X);
            if (isempty(p))
                output = ones(size(X));
            else
                output = reshape(prod(double(gammaz(pp+xx.*XX)), 2), size(X));
            end
        elseif(nargin==5)
            [~, XX] = meshgrid(p, X);
            xx = meshgrid(x,X);
            yy = meshgrid(y, X);
            [pp, YY] = meshgrid(p, Y);
            if (isempty(p))
                output = ones(size(X));
            else
                output = reshape(prod(double(gammaz(pp+xx.*XX+yy.*YY)), 2), size(X));
            end
        else
            error('nargin');
        end
    end

end
