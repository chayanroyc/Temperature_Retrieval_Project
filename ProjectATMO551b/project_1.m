% Project 1 (Temperature Retrieval)

clear all;

% initialize
nz= 51; % number of vertical layers
dz= 1;  % layer thickness
z = [0:dz:nz-1]; % vector of vertical layers (height) from 0 to 51 km

% get true x
load x_rad.dat;
xt = x_rad; 
nx = length(xt);

% get radiance obs
load y_rad.dat;
y=y_rad;
ny = length(y);

% get weighting function for each channel
load h1_rad.dat;
load h2_rad.dat;
load h3_rad.dat;
load h4_rad.dat;
load h5_rad.dat;
load h6_rad.dat;
load h7_rad.dat;
load h8_rad.dat;
load h9_rad.dat;
load h10_rad.dat;

% put them in a matrix (row by row)
H=[h1_rad;h2_rad;h3_rad;h4_rad;h5_rad; ...
   h6_rad;h7_rad;h8_rad;h9_rad;h10_rad];
[nrows, ncols]=size(H);

%==========================================================================
% plot H first and make sense of it
figure(1)
clf(1)

% line colors 
line=[0.7, 0.8, 1.0; ...
      0.5, 0.5, 0.5; ...
      1.0, 0.6, 0.8; ...
      0.0, 0.75, 0.75; ...
      0.75, 0.0, 0.75; ...
      0.75, 0.75, 0.0; ...
      0.25, 0.25, 0.25; ...
      0.0, 0.0, 1.0; ...
      0.0, 0.5, 0.0; ...
      1.0, 0.0, 0.0];

% plot weighting function for each channel
for i=1:nrows
    h2=plot(squeeze(H(i,:))',z,'LineWidth',3);
    set(h2,'Color',line(i,:));
    hold on;
end
set(gca,'Fontsize',18);
xlabel('Weighting Function','Fontsize',18);
ylabel('Height (km)','Fontsize',18);
set(gca,'TickLength',[0.05 0.05]);
legend toggle;
hl=legend({'1','2','3','4','5','6','7','8','9','10'});
grid on;
%==========================================================================

% now let's do some retrieval calculations (estimate a temperature profile
% given radiances data for 10 channels

% get guess temperature profiles (for problem 3) this is our prior profile
load xb_rad1.dat;
load xb_rad2.dat;


% main loop (we have 10 cases) 
ncases=10; % 1: problem 1, 2: problem 3a, 3: problem 3b, 4: problem 4a
           % 5: problem 4b, 6: problem 5a, 7: problem 5b, 8: problem 5c
           % 9: method 1, 10: method 2
          
for icase=1:ncases
    
    % set up R (error covariance for y)
    if (icase==4)
        eo = 10;                        % for problem 4: a case where error in y is larger
        R = diag ( ((eo/100).*y).^2);   % diagonal error covariance in y (uncorrelated channels)
        save R_prob4a.dat R -ASCII;

    elseif (icase==5)
        eo = 20;                        % for problem 4: a case where error in y is even larger
        R = diag ( ((eo/100).*y).^2);   % diagonal error covariance in y (uncorrelated channels)
        save R_prob4b.dat R -ASCII;
    else
        eo = 5;                         % percent relative error in y (problem 1)
        R = diag ( ((eo/100).*y).^2);   % diagonal error covariance in y (uncorrelated channels)
        save R_default.dat R -ASCII;
    end

    % set up P (error covariance for x)
    Pb = zeros(nz,nz);
    if (icase==6) 
        eb = 1;         % problem 5: a case where error in x is reall small
        L=3;
        for i=1:nz
            Pb(i,i) = eb^2;
            for j=1:i
                Pb(i,j)=sqrt(Pb(i,i)*Pb(j,j))*exp(- (z(i)-z(j))^2/L^2); % apply equation 
                Pb(j,i)=Pb(i,j);                                        % since it's symmetric 
            end
        end
        save Pb_prob5a.dat Pb -ASCII;

    elseif (icase==7)
        eb = 20;        % problem 5: a case where error in x is smaller
        L=3;
        for i=1:nz
            Pb(i,i) = eb^2;
            for j=1:i
                Pb(i,j)=sqrt(Pb(i,i)*Pb(j,j))*exp(- (z(i)-z(j))^2/L^2); % apply equation 
                Pb(j,i)=Pb(i,j);                                        % since it's symmetric 
            end
        end
        save Pb_prob5b.dat Pb -ASCII;
    else
        eb = 50;        % percent relative error in y (problem 1)
        if (icase==8)
            L=5;        % problem 5: a case where length scale is longer
            for i=1:nz
                Pb(i,i) = eb^2;
                for j=1:i
                    Pb(i,j)=sqrt(Pb(i,i)*Pb(j,j))*exp(- (z(i)-z(j))^2/L^2); % apply equation 
                    Pb(j,i)=Pb(i,j);                                        % since it's symmetric 
                end
            end
            save Pb_prob5c.dat Pb -ASCII;            
        else
            L=3;        % error correlation length scale in km
            for i=1:nz
                Pb(i,i) = eb^2;
                for j=1:i
                    Pb(i,j)=sqrt(Pb(i,i)*Pb(j,j))*exp(- (z(i)-z(j))^2/L^2); % apply equation 
                    Pb(j,i)=Pb(i,j);                                        % since it's symmetric 
                end
            end
            save Pb_default.dat Pb -ASCII;
        end

    end % icases
    

    % set up xb (guess profile)
    if (icase==2)
        xb = xb_rad1;   % problem 3: a case when your guess profile is different
    elseif (icase==3)
        xb = xb_rad2;   % problem 3: a case when your guess profile is different
    else
        xb = ones(nz,1)*260; % constant guess/prior profile of 260 K (for problem 1)
    end

    if (icase<=8)
        % find K of the solution
        K = Pb*H'*inv(R+H*Pb*H');

        % calculate averaging kernel A
        A = K*H;

        % estimate or retrieval
        xa = xb + K*(y-H*xb); 

        % estimate of error (not necessary for this project) 
        Pa = Pb - K*H*Pb;      
    elseif (icase==9)
        xa = inv(H'*H)*H'*y;
    elseif (icase==10)
        xa = inv(H'*H)*H'*inv(R)*y;
    end

    %======================================================================
    % now plot results

    h=figure(10+icase);
    clf(10+icase);
    set(h,'Position',[72 120 963 618]);
    subplot(1,2,1);
    plot(xb,z,'r-','Linewidth',2);
    hold on;
    plot(xa,z,'b-','Linewidth',2);
    hold on;
    plot(xt,z,'k-','Linewidth',2);
    hl=legend({'prior','retrieval','true'});
    set(gca,'Fontsize',18);
    xlabel('Temperature (K)','Fontsize',18);
    ylabel('Height (km)','Fontsize',18);
    xlim([180 300]);
    ylim([0 50]);
    
    set(gca,'TickLength',[0.05 0.05]);
    legend toggle;
    hl=legend({'Guess','Retrieval','True'});
    grid on;
    title(['Temperature Profiles for Case: ',num2str(icase)]); 


    subplot(1,2,2);
    plot(xb-xt,z,'r-','Linewidth',2);
    hold on;
    plot(xa-xt,z,'b-','Linewidth',2);
    hold on;
    hl=legend({'prior','retrieval'});
    set(gca,'Fontsize',18);
    xlabel('Temperature Bias (K)','Fontsize',18);
    ylabel('Height (km)','Fontsize',18);
    xlim([-50 50]);
    ylim([0 50]);
    
    set(gca,'TickLength',[0.05 0.05]);
    legend toggle;
    hl=legend({'Guess','Retrieval'});
    grid on;
    title(['Temperature Bias for Case: ',num2str(icase)]);     

    if (icase<=8)
        figure(20+icase);
        clf(20+icase);

        for i=1:nrows
            h2=plot(squeeze(A(i,:))',z,'LineWidth',3);
            set(h2,'Color',line(i,:));
            hold on;
        end
        set(gca,'Fontsize',18);
        xlabel('Averaging Kernel (unitless)','Fontsize',18);
        ylabel('Height (km)','Fontsize',18);
        xlim([-0.2 0.2]);
        ylim([0 50]);

        set(gca,'TickLength',[0.05 0.05]);
        legend toggle;
        hl=legend({'1','2','3','4','5','6','7','8','9','10'});
        grid on;
        title(['Averaging Kernel for Case: ',num2str(icase)]); 
    end
    
end % cases