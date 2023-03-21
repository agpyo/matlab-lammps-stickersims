clear; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spec_str = {'6A2B2'};
spec_str2 = {'6A2B2'};
T_v = {210:10:250};
Rep = 5; fpath = 'L24/'; xLen = 200; Abox = 30*30;
M = 125; thr = 0.5; thrA = 0.2; gs = 0.1;
Pf = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = length(spec_str); Mq = round(M/8);
edref = linspace(-xLen/2,xLen/2,M+1);
xref = edref(1:end-1) + 0.5*(edref(2)-edref(1));
Wout = cell(L,1); Wout2 = cell(L,1); rho = cell(L,1);
for jj = 1:L
    LT = length(T_v{jj}); sp = spec_str{jj}; sp2 = spec_str2{jj};
    fprintf(['Starting ', sp,' ...\n']);
    WT = zeros(LT,3); Wout2T = zeros(LT,5);
    rhoT = zeros(LT,2);
    for kk = 1:LT       
       T = T_v{jj}(kk);
       Wrep = zeros(Rep,1); Wtemp = zeros(1e3*Rep,2);
       rhoR = zeros(Rep,1);
       for rr = 1:Rep
           cc = 0; yall = zeros(1e2,1);
           fn = ['XYZFILES_', sp, '_run2/L_24_N_625_', sp2, '_T' , ...
                num2str(T),'_Rep', num2str(rr) ,'_POS.mat'];
           load([fpath,fn],'XYZ'); R = length(XYZ(1,1,:));
           for qq = 1:R
               A = XYZ(:,2, qq); yh0 = histcounts(A, edref); yh = yh0';
               yhs = smooth(yh0,'loess'); ys = sort(yhs); ymax = mean(ys(end-Mq:end));
               imax1_2 = find(yh > thr*ymax,1); imax1_1 = imax1_2 - 1;
               imax2_1 = find(yh > thr*ymax,1,'last'); imax2_2 = imax2_1 + 1;
               if imax1_1 > 0.125*M && imax2_2 < 0.875*M
                   x1 = xref(imax1_1) + (thr*ymax-yh(imax1_1))...
                       *(xref(imax1_2)-xref(imax1_1))/(yh(imax1_2)-yh(imax1_1));
                   x2 = xref(imax2_1) + (thr*ymax-yh(imax2_1))...
                       *(xref(imax2_2)-xref(imax2_1))/(yh(imax2_2)-yh(imax2_1));
                   imid1 = find(abs(xref-(0.5*(x1+x2)+gs*(x2-x1))) == min(abs(xref - (0.5*(x1+x2)+gs*(x2-x1)))),1);
                   imid2 = find(abs(xref-(0.5*(x1+x2)-gs*(x2-x1))) == min(abs(xref - (0.5*(x1+x2)-gs*(x2-x1)))),1);
                   dX = x2-x1; iL = find(xref < x1-dX/1.75,1, 'last');
                   iR = find(xref > x2+dX/1.75,1);
                   if ~isempty(iL) && ~isempty(iR)
                       yh1 = yh(iL:imid1); yh2 = yh(imid2:iR);
                       xh1 = xref(iL:imid1) - x1; xh2 = -(xref(imid2:iR) - x2);
                       ymin = 0.5*(mean(yh(1:iL))+mean(yh(iR:end)));
                       yf1 = 2*(yh1 - ymin)/(ymax-ymin)-1;
                       yf2 = 2*(yh2 - ymin)/(ymax-ymin)-1;
                       [f1, gof1] = cFT(xh1,yf1); [f2, gof2] = cFT(xh2,yf2);
                       if Pf == 1 && any(qq == 1:25:R)
                           xff1 = linspace(xh1(1),xh1(end),1e3)';
                           xff2 = linspace(xh2(1),xh2(end),1e3)';
                           figure(2); hold off;subplot(121);
                           plot(xh1,yf1,'k.-'); hold on; plot(xh2,yf2,'r.-');
                           plot(xff1, f1(xff1),'k'); plot(xff2, f2(xff2),'r');
                           xlim([xh1(1),xh1(end)]); ylim([-1.5,1.5])
                           drawnow
                       end
                       W1 = f1.b; W2 = f2.b; err = abs(W1-W2)/min([W1,W2]);
                       if err < thrA
                           cc = cc + 1;
                           Wtemp(cc,:) = [W1, W2];
                           if cc == 1
                               xfit = linspace(xh1(1),xh1(end),1e2)'; 
                           end
                           fy1 = fit(xh1',yh1,'linear'); fy2 = fit(xh2',yh2,'linear');
                           yall = yall + fy1(xfit) + fy2(xfit);
                          
                       end
                   end
               end
           end
           Wp = rmoutliers(reshape(Wtemp(1:cc,:),1,[])');
           WT(kk,:) = [mean(Wp), std(Wp), length(Wp)];
           
           rL = mean(yall(1:15)); rR =mean(yall(85:end));
           rhoR(rr) = (rR - rL)/(Abox*cc);
           yall2 = 2*(yall-rL)/(rR-rL)-1;
           [frs, gof] = cFT(xfit',yall2);
           if Pf == 1
              figure(2); hold off; subplot(122); plot(xfit,yall2,'.k');
              hold on; plot(xfit,frs(xfit),'r'); 
              ylim([-1.5,1.5]); xlim([xfit(1),xfit(end)]);
              drawnow;
           end 
           Wrep(rr) = frs.b;
       end
       Wout2T(kk,:) = Wrep';
       rhoT(kk,:) = [mean(rhoR),std(rhoR)/sqrt(Rep-1)];
    end
    Wout{jj} = WT; Wout2{jj} = Wout2T; rho{jj} = rhoT;
end


function [f, gof] = cFT(x,y)
% Set up fittype and options.
%ft = fittype( 'tanh(2*(x-a)/b)', 'independent', 'x', 'dependent', 'y' );
ft = fittype( 'sqrt(2)*tanh((x-a)/(2*b))/(3-tanh((x-a)/(2*b))^2)^(1/2)', 'independent', 'x', 'dependent', 'y' );
%ft = fittype( 'tanh((x-a)/(2*b))*(1-(pi/(6*sqrt(3)))*sech((x-a)/(2*b))^2)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.startpoint = [0, 5];
% Fit model to data.
[f, gof] = fit(x',y, ft, opts );
end