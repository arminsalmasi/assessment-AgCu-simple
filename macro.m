%L0=a0+b*0T
%L1=a1+b1*T
%g=(1-x)*G0Ag+x*G0+R*T*((1-x)*ln(1-x)+xln(x))+(x-x^2)*(L0+(1-2*x)*L1);
%Dg=G0Cu-G0Ag+R*T*(log(xCu)-log(1-xCu))+L0+L1-xCu*(2*L0+6*L1)+6*L1*xCu^2;
%DDG=R*T*((1/1-xCu)+(1/(xCu)))+12*L1*xCu-(2*L0+6*L1);

close all
clear variables
clc

load('activity.txt')
load('exp.txt')

syms a0FCC b0FCC a1FCC b1FCC a0Liq b0Liq a1Liq b1Liq


R = 8.314;
%% Standard Gibbs free energies
% given in the problem as an assumption
G0AgFCC = 0 ;
G0CuFCC = 0;
% from heat of melting, knowing that Sm=Hm/Tm and GL = Hm - TSm
TmAg = 1235;       % melting point of Ag
TmCu = 1358;       % melting point of Cu

HmAg = 11300;      % enthalpy of melting is J/mol for silver
HmCu = 13300;      % enthalpy of melting is J/mol for copper

%% Eutactic point
Teu = 1053; % Eutectic t
xCuEuLiq = 0.4; % Eutactic composition
xCuFCCEu1 = 0.13;
xCuFCCEu2 = 0.955;
xCuLiqEu  = 0.4;
DgFCCEu1 = G0CuFCC-G0AgFCC + ...
           +R*Teu*(log(xCuFCCEu1)-log(1-xCuFCCEu1)) + ...
           (a0FCC+Teu*b0FCC)+(a1FCC+Teu*b1FCC) - ...
           xCuFCCEu1*(2*(a0FCC+Teu*b0FCC)+6*(a1FCC+Teu*b1FCC)) +...
           6*(a1FCC+Teu*b1FCC)*xCuFCCEu1^2 ;
    
DgFCCEu2 = G0CuFCC-G0AgFCC + ...
           +R*Teu*(log(xCuFCCEu2)-log(1-xCuFCCEu2)) + ...
           (a0FCC+Teu*b0FCC)+(a1FCC+Teu*b1FCC) - ...
           xCuFCCEu2*(2*(a0FCC+Teu*b0FCC)+6*(a1FCC+Teu*b1FCC)) +...
           6*(a1FCC+Teu*b1FCC)*xCuFCCEu2^2;
        
G0AgLiq = HmAg - Teu * HmAg/TmAg;
G0CuLiq = HmCu - Teu * HmCu/TmCu;
DgLiqEu = G0CuLiq-G0AgLiq + ...
        +R*Teu*(log(xCuLiqEu)-log(1-xCuLiqEu)) + ...
        (a0Liq+Teu*b0Liq)+(a1Liq+Teu*b1Liq) - ...
        xCuLiqEu*(2*(a0Liq+Teu*b0Liq)+6*(a1Liq+Teu*b1Liq))+...
        6*(a1Liq+Teu*b1Liq)*xCuLiqEu^2 ;

    
gFCCEu1 = (1-xCuFCCEu1)*G0AgFCC + xCuFCCEu1*G0CuFCC+ ...
   Teu*((1-xCuFCCEu1).*log(1-xCuFCCEu1)+xCuFCCEu1.*log(xCuFCCEu1)) + ...
   (xCuFCCEu1-xCuFCCEu1.^2).*((a0FCC+Teu*b0FCC)+...
   (1-2*xCuFCCEu1)*(a1FCC+Teu*b1FCC));
gFCCEu2 = (1-xCuFCCEu2)*G0AgFCC + xCuFCCEu2*G0CuFCC+ ...
   R*Teu*((1-xCuFCCEu2).*log(1-xCuFCCEu2)+xCuFCCEu2.*log(xCuFCCEu2)) + ...
   (xCuFCCEu2-xCuFCCEu2.^2).*((a0FCC+Teu*b0FCC)+...
   (1-2*xCuFCCEu2)*(a1FCC+Teu*b1FCC));
gLiqEu = (1-xCuLiqEu)*G0AgFCC + xCuLiqEu*G0CuFCC+ ...
    R*Teu*((1-xCuLiqEu).*log(1-xCuLiqEu)+xCuLiqEu.*log(xCuLiqEu)) + ...
    (xCuFCCEu2-xCuLiqEu.^2).*((a0Liq+Teu*b0Liq)+...
    (1-2*xCuLiqEu)*(a1Liq+Teu*b1Liq));

mTangLine = (gFCCEu2-gFCCEu1)/(xCuFCCEu2-xCuFCCEu1);

firstline = (gLiqEu - mTangLine * xCuLiqEu) - ...
    (gFCCEu2 - mTangLine * xCuFCCEu2);
secondline = (gLiqEu - mTangLine * xCuLiqEu) - ...
    (gFCCEu1 - mTangLine * xCuFCCEu1);
thirdline=  (gFCCEu2 - mTangLine * xCuFCCEu2)- ...
    (gFCCEu1 - mTangLine * xCuFCCEu1);     

%% Top of miscibility gap
TMG = 1388;   % T top of the miscibility gap
xCuFCCMG = 0.63; % xCu top of the miscibility gap

DgFCCMG = G0CuFCC-G0AgFCC + ...
          +R*TMG*(log(xCuFCCMG)-log(1-xCuFCCMG)) + ...
          (a0FCC+TMG*b0FCC)+(a1FCC+TMG*b1FCC) - ...
          xCuFCCMG*(2*(a0FCC+TMG*b0FCC)+6*(a1FCC+TMG*b1FCC)) +...
          6*(a1FCC+TMG*b1FCC)*xCuFCCMG^2 ;

DDgFCCMG = R*TMG*((1/1-xCuFCCMG)+(1/(xCuFCCMG)))+...
           12*(a1FCC+TMG*b1FCC)*xCuFCCMG - ...
           (2*(a0FCC+TMG*b0FCC)+6*(a1FCC+TMG*b1FCC));

%% Activity of components in liquid
%  GEAgLiq = R*T*Ln(aAg/xAg) = R*T*Ln(GammaAg)
%  GECuLiq = R*T*Ln(aCu/xCu) = R*T*Ln(GammaCu)
%  GEmLiq = xAg*GEAgLiq+xCu*GECuLiq = (1-xCu)*GEAgLiq+xCu*GECuLiq
%  GEmLiq = xAg*xCu*(L0Liq+(xAg-xCu)*L1Liq) = 
%                         (xCu-xCu^2)*(L0Liq+(1-2*xCu)*L1Liq) =
%                         x(L1+L0)-x^2*(3*L1+L0)+2*L1*x^3 =
%                         L0*(x-x^2)+L1*(x-3*x^2+2*x^3)
TAC = 1400; %temperature in whihc activity is given
GEAgLiqAC = R* TAC * log(activity(:,4));
GECuLiqAC = R* TAC * log(activity(:,6));
GEmLiqAC = GEAgLiqAC.*activity(:,1)+ GECuLiqAC.*activity(:,2);

g=fit(activity(:,1),GEmLiqAC,'poly3'); 
x = linspace(0,1); 
cp =g.p1.*x.^3 + g.p2.*x.^2 + g.p3.*x+g.p4;

%plot GEM AND FIT
figure
box on
hold on
    h= plot(x,cp,activity(:,2),GEmLiqAC,'o');
    legendCell(1) = cellstr('$Fit \:^{Liq}G^E_m$');
    legendCell(2) = cellstr('$^{Liq}G^E_m$');
    title('$^{Liq}G^E_m$ vs $x_{Cu}$ at 1400 K','interpreter','latex');
    legend(legendCell(:),'interpreter','latex');
    ax = gca;
    ax.FontSize = 15;
    set(h,{'markers'},{2;10})  
    set(h,'LineWidth',2);
    ylim([0,3500]);
    D ='$Fit \:^{Liq}G^E_m \approx 0x^3-1.3\times10^4x^2+1.3\times10^4x+0 $'; 
    h=annotation(gcf,'textbox',[.25 .2 .2581 .0929],...
        'string',D,'interpreter','latex');
    h.FitBoxToText = 'on';
hold off

%plot activity
%order = ['xAg '	'xCu '	'aAg '	'gammaAg '	'acu '	'gammaCu '];
   
figure
box on   
hold on
    plot( activity(:,1), activity(:,3), 'b' )
        legendCell(1) = cellstr('$a_{Ag}$');
    plot( 1-activity(:,2), activity(:,5), 'g' )
        legendCell(2) = cellstr('$a_{Cu}$');
    plot( linspace(0,1),linspace(0,1), 'r-.', linspace(1,0),...
        1-linspace(1,0), 'm--' )
        legendCell(3) = cellstr('$a_{ideal}^{Ag}$');
        legendCell(4) = cellstr('$a_{ideal}^{Cu}$');
    plot( exp(:,1),exp(:,3), 'O')
        legendCell(5) = cellstr('$a_{Ag}^{Experiment}$');
    title('Activty of Ag and Cu in Liquid at 1400 K',...
        'interpreter','latex')
    ylabel('ACTIVITY','interpreter','latex');
    xlabel('$X_{Ag}$','interpreter','latex');
    legend(legendCell(:),'interpreter','latex')

    ax = gca;
    ax.FontSize = 15;
hold off

%% Solver
eq1 = DgFCCEu1-DgFCCEu2 == 0       %eutactic tangant slope
eq2 = DgFCCEu2 - DgLiqEu == 0      %eutactic tangant slope
eq3 = DgLiqEu - DgFCCEu1 == 0      %eutactic tangant slope
eq4 = DgFCCMG == 0                 %miscibility gap   
eq5 = DDgFCCMG == 0                %misibility gap
eq6 = (a1Liq+TAC*b1Liq) == 0       %activity
eq7 = (a0Liq+TAC*b0Liq) - g.p2 ==0 %activity
eq8 = firstline == 0 %eutactic tangant distance from origin
eq9 = secondline== 0%eutactic tangant distance from origin
eq10=thirdline==0 %eutactic tangant distance from origin
[A,B]=equationsToMatrix([eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9,eq10],...
                    [a0FCC,b0FCC,a1FCC,b1FCC,a0Liq,b0Liq,a1Liq,b1Liq]);
A =double(A)
B =double(B)

C2 = A\B
   
printStr= ([
' %s %1.3e%s%1.3e%s \n %s',...
' %1.3e%s%1.3e%s \n %s %1.3e%s%1.3e%s \n %s %1.3e%s%1.3e%s \n']);
fprintf( printStr,...
         'L0S =' ,C2(1,1),'+',C2(2,1), '*T',...
         'L1S = ',C2(3,1),'+', C2(4,1),'*T',...
         'L0L = ',C2(5,1),'+', C2(6,1),'*T',...
         'L1L = ',C2(7,1),'+', C2(8,1),'*T');

fprintf('%s%f%s%f \n', 'G0AgLiq = ', HmAg, '- T*', HmAg/TmAg);

fprintf('%s%f%s%f \n',  'G0CuLiq =', HmCu, '- T*', HmCu/TmCu);

%% TEST PLOT     

for T = [ 1001]%1:100:1400;
    figure
    hold on
    plotFCC=gplot(T, G0AgFCC, G0CuFCC, C2(1,1), C2(2,1), C2(3,1), C2(4,1));
    legendCell(1) = cellstr('$FCC$');
    plotLiq=gplot(T, G0AgFCC+HmAg-T*HmAg/TmAg, G0CuFCC+HmCu-T*HmCu/TmCu,...
              C2(5,1), C2(6,1), C2(7,1), C2(8,1));
    legendCell(2) = cellstr('$Liquid$');
    legendCell(5) = cellstr('$a_{Ag}^{Experiment}$');
    title('Gibbs free energies at eutectic temperature', 'interpreter','latex')
    ylabel('G','interpreter','latex');
    xlabel('$X_{Cu}$','interpreter','latex');
    box on
    legend(legendCell(:),'interpreter','latex')
  
    ax = gca;
    ax.FontSize = 15;      
          
    hold off
end
 printStr= ...
' %s  %e \n %s %e \n %s %e \n %s %e \n %s %e \n %s %e \n %s %e \n %s %e \n';
fprintf( printStr,...
         'a0FCC =' ,C2(1,1), 'b0FCC = ',C2(2,1),...
         'a1FCC = ',C2(3,1), 'b1FCC = ',C2(4,1),...
         'a0Liq = ',C2(5,1), 'b0Liq = ',C2(6,1),...
         'a1Liq = ',C2(7,1), 'b1Liq = ',C2(8,1));


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%**************************************************************************
function [a]=gplot(T, G01, G02, a0, b0, a1, b1)
    R=8.314;
    x= linspace(0,1,100);
    g = (1-x)*G01 + x*G02+ ...
          R*T*((1-x).*log(1-x)+x.*log(x)) + ...
          (x-x.^2).*((a0+b0*T)+(1-2*x)*(a1+b1*T));
    plot(x,g,'LineWidth',2)
    a=0;
end