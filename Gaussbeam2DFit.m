function Gaussbeam2DFit

pixelsize = 3.75;                                                          %Pixelsize of the camera, start parameter = 3.75





offsetOnOff = 0;distPrevIm = 0;wR = 0;

                                                                           %read in the first image
picture = readInImage;

E = createGaussFitOfImage(picture, offsetOnOff,wR);                           %create first gauss fit, return all relevant parameters in E

saveGaussFitCell {1,1} = {E{1,:}, picture, distPrevIm};                    %save all parameters in Cell array "saveGaussFitCell"



plotImFit(saveGaussFitCell {1,1});                                                      %plot first image and corresponding fit




            uicontrol('Style', 'pushbutton', 'String', 'Intensity Profile',...
                'Position', [450 120 100 20],...
                'Callback', @pushbuttonIP_Callback); 
                                                                        %pushbutton for adding more images
            uicontrol('Style', 'pushbutton', 'String', 'Add Image',...
                'Position', [450 30 100 20],...
                'Callback', @pushbutton1_Callback);                         %save data at txt
            
            uicontrol('Style', 'pushbutton', 'String', 'Save Results',...
            'Position', [450 60 100 20],...
            'Callback', @pushbutton2_Callback);  
                                                                            %checkbox for offset
         if offsetOnOff == 0,
            OffsetOnOffcheck = uicontrol('Style', 'checkbox', 'String', 'Offset',...
            'Position', [450 90 100 20],...
            'Callback', @checkbox1_Callback);
         else
         end
            
            uicontrol('Style', 'text',...
           'String', 'Pixelsize [µm]',...
           'Position', [360 60 80 20],...
           'Callback', @setmap); 
        
            uicontrol('Style', 'popup',...
           'String', {'3.75','3.10','3.69','4.54','5.55','5.86','6.45'},...
           'Position', [360 30 80 20],...
           'Callback', @popup1_callback); 
function tempImg = readInImage  
         
          
          imgfile=imgetfile;          
             tempImg = imread(imgfile);    
             imshow(tempImg, 'InitialMagnification', 50);
             uiwait(msgbox('Select a square. Rectangle will be cut off to square','Read in image','modal'));
             tempImg = imcrop();
             tempImg = im2double(tempImg);          
             close all;
                
end 
function gaussFitResult = createGaussFitOfImage(picture, offsetOnOff,wR)
                           
                           %make picture square
    if length(picture(:,1)) <=length(picture(1,:)),    
         x = (1:length(picture(:,1)));
         y = (1:length(picture(:,1)));
    else 
         x = (1:length(picture(1,:)));
         y = (1:length(picture(1,:)));
    end
    
                        %generating grid
[x,y] = meshgrid(x,y);

                        %cut picture to length of x,y
picpoint = picture(1:length(x),1:length(y));

                       
                        %find start parameters

params = findStartParameters(picpoint,wR);
      
       




        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        %fit with function / minimize
                        %set optimset max iterationsm maxfunEvals to 10000000
                        %decide wether Offset is on or off
if offsetOnOff == 0,                   
        parameters = [params(6), params(4), params(5),params(3),params(1),params(2)];
else
        parameters = [params(6), params(4), params(5),params(3),params(1),params(2),params(7)];
end

        opt = optimset ( 'MaxIter' , 100000000, 'MaxFunEvals', 100000000 );
chifunc = fminsearch(@fit2d,parameters, opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

                                                                            %use minimized data to calculat f(x,y), normalize it
                                                                            %set angel +90° if waist x < waist y
gaussFitResult = preData(chifunc);



    function gaussFitResult = preData(chifunc)


        w = chifunc(1); wbx = chifunc(2); wby = chifunc(3); Imax = chifunc(4);
        cx = chifunc(5); cy = chifunc(6); 
        
        if offsetOnOff == 0,
            offst = 0;
                                   %normalize chi square                        
            f = Imax*exp(-(2/(wbx^2))*(cos(w)*(x-cx)+sin(w)*(y-cy)).^2-(2/(wby^2))*(-sin(w)*(x-cx)+cos(w)*(y-cy)).^2);

        else
                 offst = chifunc(7);
                                  %normalize chi square                        
                f = Imax*exp(-(2/(wbx^2))*(cos(w)*(x-cx)+sin(w)*(y-cy)).^2-(2/(wby^2))*(-sin(w)*(x-cx)+cos(w)*(y-cy)).^2)-offst;
        end

    
        chi = sum(sum((f-picpoint).^2))/(sum(sum((f).^2)));
                        

                        %angle in degrees
        w = 180*w/(pi);
        wR = w;
                        %set bigger waist in x direction
                   
    %    if wbx <= wby 
     %       switchyandx = wbx;
      %      wbx = wby;
       %     wby = switchyandx;
        %    if w <+0
         %   w = w + 90;
          %  else 
           % w = w - 90;
            %end
        
  %     else    
   %     end
       
        if w >= 0
            while w>=180
                w = w-360;
            end
        else 
            while w<=-180
                w = w+360;
            end
        end

    gaussFitResult = {chi, x, y, f, wbx, wby, Imax, picpoint,offst,w,cx,cy,wR};
    end
    function chi = fit2d(parameters)                                        %Minimized function
                                                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        w = parameters(1); wbx = parameters(2); wby = parameters(3);        %width (wbx,wby) = 4sigma ====D4sigma/second moment width
        Imax = parameters(4); cx = parameters(5); cy = parameters(6);
       if offsetOnOff == 0,
           f = Imax*exp(-(2/(wbx^2))*(cos(w)*(x-cx)+sin(w)*(y-cy)).^2-(2/(wby^2))*(-sin(w)*(x-cx)+cos(w)*(y-cy)).^2);
        
        else    offst = parameters(7);
            f = Imax*exp(-(2/(wbx^2))*(cos(w)*(x-cx)+sin(w)*(y-cy)).^2-(2/(wby^2))*(-sin(w)*(x-cx)+cos(w)*(y-cy)).^2)-offst;
        
        end
        
       chi = sum(sum((f-picpoint).^2)); 
       
       
    end
    function params = findStartParameters(picpoint,wR)
        
        %sum over two dimensions of image
        profilex = sum(picpoint,1);
        profiley = sum(picpoint,2);
        
        [~,cx] = max(profilex);
        [~,cy] = max(profiley);
        Imax = picpoint(cx,cy);
        
        %offset
        offst = picpoint(1,1);
        
        %find maximum and minimum FWHM dependency on angle
        imFWHM = (picpoint - Imax/2);
        i = cx;
        j = cy;   
       
            while imFWHM(i,j) > 0&& abs(i)<length(imFWHM(1,:))&& abs(i)<length(imFWHM(:,1))&&i>1&&j>1,
                i = i+1;                        
            end   
                i1 = i;
                i = cx;
            while imFWHM(i,j) >0&& abs(i)<length(imFWHM(1,:))&& abs(i)<length(imFWHM(:,1))&&i>1&&j>1,
                i = i-1;
            end
                widthx = abs(i) + abs(i1);
                i = cx;
            while imFWHM(i,j) > 0&& abs(j)<length(imFWHM(1,:))&& abs(j)<length(imFWHM(:,1))&&i>1&&j>1,
                j = j+1;                        
            end   
                j1 = i;
                j = cy;
            while imFWHM(i,j) >0&& abs(j)<length(imFWHM(1,:))&& abs(j)<length(imFWHM(:,1))&&i>1&&j>1,
                j = j-1;
            end
                widthy = abs(j) + abs(j1);
        
        w = wR;
        
        
        params = [cx,cy,Imax, widthx, widthy, w, offst];
    end  
end    
function plotImFit(E)
        
        chi = E{1,1};x = E{1,2}; y = E{1,3}; f = E{1,4}; wbx = E{1,5};wby = E{1,6};Imax = E{1,7}; picpoint = E{1,8};
        offst = E{1,9}; w = E{1,10};
        
       
        %plot fit, with pixelsize
        subplot(222),surf(x*pixelsize,y*pixelsize,f),shading flat
            xlabel('x'); ylabel('y'); title('fit')
            pause(.01)
            axis([0,length(x)*pixelsize,0, length(x)*pixelsize, 0 Imax])
        
        %plot image
        subplot(221),surf(x*pixelsize,y*pixelsize,picpoint),shading flat
            xlabel('x'); ylabel('y'); title('image')
            axis([0, length(x)*pixelsize,0, length(x)*pixelsize, 0 Imax])
            
        %plot data table
        
            rnames = {'waist x [µm]','waist y [µm]','angle(x) [°]','maximum intensity','normalized chi-square','offset' };
            d = {(wbx*pixelsize);(wby*pixelsize);w';Imax;chi;offst};
            
            uitable('Data', d, 'RowName',rnames,'ColumnName','values','ColumnWidth','auto', 'Position', [20 30 330 140]);
        
    
end
function saveGaussFitCell = addMoreImages(saveGaussFitCell,counter)
    close all;
        counter = counter +1;
        picture = readInImage;
        E = createGaussFitOfImage(picture, offsetOnOff,wR);
        
        distance2first = inputdlg('Enter Distance to starting point [mm]:',...
             'Distance', [1 50]);
            distPrevIm = str2double(distance2first{:})*1000;
        saveGaussFitCell {counter,1}  = {E{1,:}, picture, distPrevIm};

    close all;
        figure('Name','Results');   
        
        plotImagesInRow(saveGaussFitCell, counter);
         
         uicontrol('Style', 'pushbutton', 'String', 'Add Next Image',...
            'units','Normalized','Position', [.01 .39 .15 .05],...
            'Callback', @pushbutton3_Callback);
         
         uicontrol('Style', 'pushbutton', 'String', 'Create Profile',...
            'units','Normalized','Position', [.33 .39 .15 .05],...
            'Callback', @pushbutton5_Callback);
        
        uicontrol('Style', 'pushbutton', 'String', 'New Measurement',...
            'units','Normalized','Position', [.17 .39 .15 .05],...
            'Callback', @pushbutton6_Callback);
        uicontrol('Style', 'pushbutton', 'String', 'Save All',...
            'units','Normalized','Position', [.49 .39 .15 .05],...
            'Callback', @pushbutton7_Callback);
         uicontrol('Style', 'text',...
           'String', 'Intensity Profile',...
           'units','Normalized','Position', [.65 .39 .15 .05],...
           'Callback', @setmap);
         uicontrol('Style', 'popup',...
           'String', {'1','2','3','4','5','6','7','8','9','10'},...
           'units','Normalized','Position', [.81 .39 .15 .05],...
           'Callback', {@popup2_callback,saveGaussFitCell}); 
       
     function pushbutton7_Callback(~, ~, ~)
                
            lengthCell = size(saveGaussFitCell);
            
            fileID = fopen('beamFitFileAll.txt','w');
            fprintf(fileID,'%12s %12s %12s %12s %12s %12s %12s \n','imagenumber','normalized chi-square','waist_x','waist_y','I_max','Offset','Angle');
         for sC = 1:lengthCell(1);
            saveA = saveGaussFitCell{sC,1};
            saveA_i = [sC;saveA{1,1};saveA{1,5}*pixelsize;saveA{1,6}*pixelsize;saveA{1,7};saveA{1,9};saveA{1,10}];
            fprintf(fileID,'%6.0f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n',saveA_i);
         end
        fclose(fileID);   
    end        
     function pushbutton5_Callback(~, ~, ~)
       
       giveToTable = createGaussProfile(saveGaussFitCell);
       
       uicontrol('Style', 'pushbutton', 'String', 'Save Results',...
            'units','Normalized','Position', [.01 .26 .2 .05],...
            'Callback', @pushbutton4_Callback);
        
        function pushbutton4_Callback(~, ~, ~)                        
            saveA = saveGaussFitCell{1,1};
            saveA_i = [1;saveA{1,1};saveA{1,5}*pixelsize;saveA{1,6}*pixelsize;saveA{1,7};saveA{1,9};saveA{1,10};giveToTable{1,1};giveToTable{1,2};giveToTable{2,2};giveToTable{1,3};giveToTable{2,3};giveToTable{1,4};giveToTable{2,4};giveToTable{1,5};giveToTable{2,5};giveToTable{1,6};giveToTable{2,6}];
         fileID = fopen('beamFitFile1.txt','w');
        fprintf(fileID,'%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n','imagenumber','normalized chi-square','waist_x','waist_y','I_max','Offset','Angle','wavelength','w_0_x','w_0_y','z_0_x','z_0_y','chi_square_px','chi_square_py','M_squared_x','M_squared_y','Rayleight_length_x','Rayleight_length_y');
        fprintf(fileID,'%6.0f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\r\n',saveA_i);
        fclose(fileID);   
     end       
     end          
     function pushbutton3_Callback(~, ~, ~)
        addMoreImages(saveGaussFitCell,counter);
     end
     function plotImagesInRow(saveGaussFitCell, counter)
        
        close all;
        figure('Name','Results');
        i_1 = 1;
        while i_1 <= counter,
            k = counter + i_1;
            A = saveGaussFitCell{i_1,1};
            
            chi = A{1,1}; x = A{1,2}; y=A{1,3}; wbx = A{1,5}; wby = A{1,6}; Imax = A{1,7};
            picpoint = A{1,8}; offst = A{1,9}; w = A{1,10}; 
            cx = A{1,11}; cy = A{1,12};picture = A{1,13}; distPrevImN = A{15}/1000;
            
            f = Imax*exp(-(2/(wbx^2))*(cos(w)*(x-cx)+sin(w)*(y-cy)).^2-(2/(wby^2))*(-sin(w)*(x-cx)+cos(w)*(y-cy)).^2)-offst;
            
            subplot(4,counter,k),surf(x*pixelsize,y*pixelsize,f),shading flat
                xlabel('x'); ylabel('y'); title('fit')
                axis([0, length(x)*pixelsize,0, length(x)*pixelsize, 0, Imax])
            
            subplot(4,counter,i_1),surf(x*pixelsize,y*pixelsize,picpoint),shading flat
                xlabel('x'); ylabel('y'); title('image')
                axis([0, length(x)*pixelsize,0, length(x)*pixelsize, 0, Imax])
                        

            rnames = {'waist x [µm]','waist y [µm]','angle(x) [°]','maximum intensity','normalized chi-square','offset', 'distance [mm]' };
            d{1,i_1} = (wbx*pixelsize);d{2,i_1} = (wby*pixelsize);d{3,i_1} = w; 
            d{4,i_1} = Imax;d{5,i_1} = chi;d{6,i_1} = offst; d{7,i_1} = distPrevImN;
            
            
            uitable('Data', d, 'RowName',rnames,'ColumnWidth','auto','ColumnName','numbered', 'units','Normalized','Position', [.01 .01 .9 .37]);
            i_1 = i_1 + 1;
        end
        
    end
end
function giveToTable = createGaussProfile(saveGaussFitCell)
        
        
        wavInput = inputdlg('Enter the wavelength [nm]:',...               %Input Wavelength
             'Distance', [1 50]);
             wavelength = str2double(wavInput{:});
             wavelength = wavelength/(1000);
             
        figure('Name','Profiles')
        o = 5;
        while o <= 6                                                       %Do profile for x and y       
            width_wbxy = 0;
            min_width = 0;
            lengthCell = size(saveGaussFitCell);
            n_i = 1;
            while n_i <= lengthCell(1)                                      %run over all images
                
                 A = saveGaussFitCell{n_i,1};
                 if n_i == 1
                     addDist = 0; z = 0;
                 else
                 addDist = A{1,15};
                 z = [z,addDist];
                 end          
                       if n_i == 1,
                           width_wbxy = [A{1,o}];
                           
                       else
                           width_wbxy = [width_wbxy, A{1,o}];               %set waist for x and y of the current image to 'width_wbxy'
                       end
                     
                     if A{1,o} <= min_width,
                         min_width = A{1,o};
                     elseif min_width == 0,
                         min_width = A{1,o};
                         else
                     end                   
                     
                n_i = n_i +1;
            end
            wZ = min_width*pixelsize;                                       % Starting parameter for fit wZ set to minimal waist of x or y images
            z_0 = width_wbxy(1,2);                                          %set starting parameter for z0 on position of image 2                                                 
            M_square = 5;
            proParameters = [wZ, z_0, M_square];
            
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %find minimum for w_0 and z_0                                                      
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
Aeq=[];
beq=[];
lb=[-inf,-inf,1];
ub=[];
options=optimset('Algorithm','active-set','MaxIter' , 1000000, 'MaxFunEvals', 1000000);

            min_wZ = fmincon(@proParams,proParameters,[],[],Aeq,beq,lb,ub,{},options);
            
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            z_xy = 0:0.01:z(length(z));            
            width_wb = min_wZ(1) * sqrt( 1+ min_wZ(3)^4.*((z_xy-min_wZ(2)) / (pi*(min_wZ(1).^2)/wavelength)).^2);
            
            wwb1 = min_wZ(1) * sqrt( 1+ min_wZ(3)^4.*((z-min_wZ(2)) / (pi*(min_wZ(1).^2)/wavelength)).^2);
            chi_square_p = sum((wwb1 - width_wbxy*pixelsize).^2)/(sum((wwb1).^2));
            k = o - 4;
            subplot(3,1,k),plot(z_xy/1000, width_wb) 
                hold on;
                plot(z/1000,width_wbxy*pixelsize,'or')
            if o == 5,
                xlabel('z [mm]'); ylabel('w(z) [µm]'); title('Profile waist in x');
            else
                xlabel('z [mm]'); ylabel('w(z) [µm]'); title('Profile waist in y');
            end
            o = o + 1; 
            
            ray_leng = pi*((min_wZ(1))^2)/wavelength/1000;
            
            giveToTable{k,1} = wavelength*1000; giveToTable{k,2} = min_wZ(1); giveToTable{k,5} = min_wZ(3);
            giveToTable{k,3} = min_wZ(2)/1000; giveToTable{k,4} = chi_square_p; giveToTable{k,6}=ray_leng;
        end
            hold off;
     
            
                      
            rnames = {'wavelength [nm]','w_0 [µm]','z_0 [mm]', 'normalized chi-square','M-square','Rayleigh length [mm]'};
            d{1,1} = giveToTable{1,1}; d{2,1} = giveToTable{1,2}; d{3,1} = giveToTable{1,3}; d{4,1} = giveToTable{1,4};d{5,1} = giveToTable{1,5};
            d{1,2} = giveToTable{2,1}; d{2,2} = giveToTable{2,2}; d{3,2} = giveToTable{2,3}; d{4,2} = giveToTable{2,4};d{5,2} = giveToTable{2,5};
            d{6,1} = giveToTable{1,6};d{6,2} = giveToTable{2,6};
            uitable('Data', d, 'RowName',rnames,'ColumnName',{'x waist', 'y waist'},'units', 'Normalized', 'Position', [.01 .01 .9 .25]);
            
            
                    
    function p = proParams(proParameters)
        
        width_wb = proParameters(1).* sqrt( 1+ proParameters(3)^4*((z-proParameters(2))./(pi*(proParameters(1)^2)/wavelength)).^2);
        p = sum((width_wb - width_wbxy*pixelsize).^2);
        
    end  
    
end
function makeProfiles(saveGaussFitCell,counterProfile)
    figure('Name','Intensity measurement and fit')
    profileVars = saveGaussFitCell{counterProfile,1};
    
    wP = profileVars{1,13};
    for fP= 1:4;
    
    xMax(fP) = profileVars{1,11};
    yMax(fP) = profileVars{1,12};
    
    W(1) = 0;
    W(2) = 0;
    W(3) =  90;
    W(4) = 90;
    
    while yMax(fP)>2 && xMax(fP)>2 && xMax(fP)<(length(profileVars{1,2})) && yMax(fP)<(length(profileVars{1,3}))
        xMax(fP) = xMax(fP) + cosd(wP+W(fP))*(-1)^(fP+sind(W(fP)));
        yMax(fP) = yMax(fP) + cosd(90-wP+W(fP))*(-1)^(fP+sind(W(fP)));
              
    end

    end
    
    for PP = 1:2;
        
    xVal = [xMax(2*PP-1) xMax(2*PP)];
    yVal = [yMax(2*PP-1) yMax(2*PP)];
    imPoint = profileVars{1,8};
    IP{PP} = improfile(imPoint,xVal, yVal);
        
    end
   
    subplot(211),
            
    fx = profileVars{1,7}*exp(-(2/profileVars{1,5}^2)*(profileVars{1,2}-profileVars{1,11}).^2)+profileVars{1,9};
    plot(profileVars{1,2},fx,'Color','blue')
    xlabel('x(w)'); ylabel('Intensity'); title('x-waist direction')
    hold on
    plot(IP{1},'Color','red')    
    hold off
    
    subplot(212),xlabel('x'); ylabel('y'); title('fit')
            
    fy = profileVars{1,7}*exp(-(2/profileVars{1,6}^2)*(profileVars{1,3}-profileVars{1,12}).^2)+profileVars{1,9};
    plot(profileVars{1,3},fy,'Color','blue')
    xlabel('y(w)'); ylabel('Intensity'); title('y-waist direction')
    hold on
    plot(IP{2},'Color','red')
    hold off
end
function pushbutton1_Callback(~, ~, ~)
        
        addMoreImages(saveGaussFitCell,1);
      
end
function pushbuttonIP_Callback(~, ~, ~)
        
        makeProfiles(saveGaussFitCell,1);
      
end
function popup2_callback(source,~,saveGaussFitCell)
                
        val1 = get(source,'Value');
        maps1 = get(source,'String'); 
        
        counterProfile1 = str2double(maps1(val1));
        
        makeProfiles(saveGaussFitCell,counterProfile1);
        
end
function pushbutton2_Callback(~, ~, ~)
    
    
        saveCounter = 1; 
        
        saveA = saveGaussFitCell{1,1};
        

            saveA_i = [saveCounter;saveA{1,1};saveA{1,5}*pixelsize;saveA{1,6}*pixelsize;saveA{1,7};saveA{1,9};saveA{1,10}];
           
            fileID = fopen('beamFitFile1.txt','w');
        fprintf(fileID,'%12s %12s %12s %12s %12s %12s %12s\n','imagenumber','normalized chi-square','waist_x','waist_y','I_max','Offset','Angle');
        fprintf(fileID,'%6.0f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\r\n',saveA_i);
        fclose(fileID);
        
     
end
function pushbutton6_Callback(~, ~, ~)
         
        clearvars;
        close all;
        Gaussbeam2DFit;
end
function checkbox1_Callback(hObject, ~, ~)
                if (get(hObject,'Value') == get(hObject,'Max')),
                    display('Offset On');
                    offsetOnOff = 1;
                    OffsetOnOffcheck = 1;
                    E = createGaussFitOfImage(picture, offsetOnOff,wR);
                    saveGaussFitCell{1,1} = {E{1,:}, picture, distPrevIm};

                    plotImFit(saveGaussFitCell {1,1});
                else
                    display('Offset Off');
                    offsetOnOff = 0;
                    OffsetOnOffcheck = 0;
                    E = createGaussFitOfImage(picture, offsetOnOff,wR);
                    saveGaussFitCell{1,1} = {E{1,:}, picture, distPrevIm};
                    plotImFit(saveGaussFitCell {1,1});
                end
end
function popup1_callback(source,~)
                
        val = get(source,'Value');
        maps = get(source,'String'); 
        pixelsize = str2double(maps(val));
        E = createGaussFitOfImage(picture, offsetOnOff,wR);
        saveGaussFitCell{1,1} = {E{1,:}, picture, distPrevIm};

        plotImFit(saveGaussFitCell {1,1});
end

end