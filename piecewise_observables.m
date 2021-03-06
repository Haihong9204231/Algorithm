%状态向量由每个截面的截面应变的一个6D矢量组成。每个横截面的位置和速度等积分量的计算公式如下：
function status = piecewise_observables(t,z,flag)
global gv

X           =gv.X;
npie        =gv.npie;
disc        =gv.disc;

%----------------------------------------------------------------------
% observables: position (g), velocity (eta)

switch flag
    case []
        for zz=1:length(t)
            g                =gv.g;
            eta              =gv.eta;
            nstep            =gv.nstep;  %%%%%
            
            Xi              =z(1:6*npie,zz);
            Xidot           =z(6*npie+1:12*npie,zz);
            g_r              =[0 -1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1]; %%    % cantilever
            %the transformation between the spatial frame and the base frame of the soft manipulator
            g_prec           =diag([1 1 1 1]); %初始位形
            eta_prec         =zeros(6,1);

            for jj=1:npie
                xin             =Xi(6*(jj-1)+1:6*(jj-1)+6,:);
                xidotn          =Xidot(6*(jj-1)+1:6*(jj-1)+6,:);
                kn               =xin(1:3);
                thetan           =sqrt(kn'*kn);
    
                % kinematics
                for ii=1:disc   %% disc数量
                    invAdjgn_here    =piecewise_invAdjoint(X(ii),thetan,xin);
                    intdAdjgn_here   =piecewise_ADJ(X(ii),thetan,xin);
                    g(4*(nstep-1)+1:4*(nstep-1)+4,4*(jj-1)*disc+4*(ii-1)+1:4*(jj-1)*disc+4*(ii-1)+4)...
                        =g_r*g_prec*piecewise_expmap(X(ii),thetan,xin);
                    eta(6*(nstep-1)+1:6*(nstep-1)+6,(jj-1)*disc+ii)...
                        =invAdjgn_here*(eta_prec+intdAdjgn_here*xidotn);
                end
    
                % recursive factors
                invAdjgn_last   =piecewise_invAdjoint(X(disc),thetan,xin);
                intdAdjgn_last  =piecewise_ADJ(X(disc),thetan,xin);
                g_prec          =g_prec*piecewise_expmap(X(disc),thetan,xin);%%位置递归
                ADxin           =intdAdjgn_last*xidotn;
                eta_prec        =invAdjgn_last*(eta_prec+ADxin);  %%速度递归
            end
            
            gv.g             =g;
            gv.eta           =eta;
            gv.nstep         =nstep;
            gv.nstep         =nstep+1;
        end
end
  
status      =0;  %%

% eof