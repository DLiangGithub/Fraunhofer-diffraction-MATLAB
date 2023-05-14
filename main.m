classdef myapp < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure           matlab.ui.Figure
        GridLayout         matlab.ui.container.GridLayout
        LeftPanel          matlab.ui.container.Panel
        TabGroup           matlab.ui.container.TabGroup
        Tab                matlab.ui.container.Tab
        nmLabel            matlab.ui.control.Label
        mLabel             matlab.ui.control.Label
        df_la_slider       matlab.ui.control.Slider
        df_a_slider        matlab.ui.control.Slider
        df_2d_button       matlab.ui.control.Button
        df_la_spin         matlab.ui.control.Spinner
        df_a_spin          matlab.ui.control.Spinner
        Label              matlab.ui.control.Label
        df_DropDown        matlab.ui.control.DropDown
        df_dh_button       matlab.ui.control.Button
        df_1d_button       matlab.ui.control.Button
        df_3d_button       matlab.ui.control.Button
        Tab_2              matlab.ui.container.Tab
        duo_dh_button      matlab.ui.control.Button
        nmLabel_2          matlab.ui.control.Label
        mLabel_2           matlab.ui.control.Label
        duo_la_slider      matlab.ui.control.Slider
        duo_a_slider       matlab.ui.control.Slider
        duo_ty_button      matlab.ui.control.Button
        duo_la_spin        matlab.ui.control.Spinner
        duo_a_spin         matlab.ui.control.Spinner
        DropDown_2Label    matlab.ui.control.Label
        duo_DropDown       matlab.ui.control.DropDown
        mmLabel            matlab.ui.control.Label
        duo_d_slider       matlab.ui.control.Slider
        duo_d_spin         matlab.ui.control.Spinner
        Label_2            matlab.ui.control.Label
        duo_n_slider       matlab.ui.control.Slider
        duo_n_spin         matlab.ui.control.Spinner
        duo_1d_button      matlab.ui.control.Button
        duo_3d_button      matlab.ui.control.Button
        Tab_3              matlab.ui.container.Tab
        nmLabel_3          matlab.ui.control.Label
        mLabel_3           matlab.ui.control.Label
        jk_la_slider       matlab.ui.control.Slider
        jk_a_slider        matlab.ui.control.Slider
        jk_2d_button       matlab.ui.control.Button
        jk_la_spin         matlab.ui.control.Spinner
        jk_a_spin          matlab.ui.control.Spinner
        jk_dh_button       matlab.ui.control.Button
        Label_3            matlab.ui.control.Label
        jk_DropDown        matlab.ui.control.DropDown
        mLabel_5           matlab.ui.control.Label
        jk_b_slider        matlab.ui.control.Slider
        jk_b_spin          matlab.ui.control.Spinner
        jk_1d_button       matlab.ui.control.Button
        jk_3d_button       matlab.ui.control.Button
        Tab_4              matlab.ui.container.Tab
        nmLabel_4          matlab.ui.control.Label
        mLabel_4           matlab.ui.control.Label
        yk_la_slider       matlab.ui.control.Slider
        yk_a_slider        matlab.ui.control.Slider
        yk_2d_button       matlab.ui.control.Button
        yk_la_spin         matlab.ui.control.Spinner
        yk_a_spin          matlab.ui.control.Spinner
        yk_dh_button       matlab.ui.control.Button
        Label_4            matlab.ui.control.Label
        yk_DropDown        matlab.ui.control.DropDown
        yk_1d_button       matlab.ui.control.Button
        yk_3d_button       matlab.ui.control.Button
        Tab_5              matlab.ui.container.Tab
        Label_5            matlab.ui.control.Label
        Label_6            matlab.ui.control.Label
        set_p_slider       matlab.ui.control.Slider
        set_a_slider       matlab.ui.control.Slider
        set_p_spin         matlab.ui.control.Spinner
        set_a_spin         matlab.ui.control.Spinner
        set_button         matlab.ui.control.Button
        set_Label          matlab.ui.control.Label
        set_defult_Button  matlab.ui.control.Button
        Label_7            matlab.ui.control.Label
        axis_switch        matlab.ui.control.Switch
        RightPanel         matlab.ui.container.Panel
        UIAxes             matlab.ui.control.UIAxes
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.UIFigure.Name = "夫琅禾费衍射仿真工具";
            global N;
            global M;
            N = 500;
            M = 20;
            axis(app.UIAxes,'on');
            % axis(app.UIAxes,'off');
        end

        % Value changed function: df_la_slider
        function df_la_sliderValueChanged(app, event)
            app.df_la_spin.Value = app.df_la_slider.Value;
        end

        % Value changed function: df_la_spin
        function df_la_spinValueChanged(app, event)
            app.df_la_slider.Value = app.df_la_spin.Value;
        end

        % Value changed function: df_a_slider
        function df_a_sliderValueChanged(app, event)
            app.df_a_spin.Value = app.df_a_slider.Value;
        end

        % Value changed function: df_a_spin
        function df_a_spinValueChanged(app, event)
            app.df_a_slider.Value = app.df_a_spin.Value;
        end

        % Button pushed function: df_1d_button
        function df_1d_buttonButtonPushed(app, event)
            global N;
            lambda = app.df_la_slider.Value*1e-9;
            a = app.df_a_slider.Value*1e-3;
            %a = (400:1:780)*1e-9;
            theta = linspace(-8e-5, 8e-5, N);  
            %theta = linspace(-5e-2, 5e-2, N);  
            x = tan(theta);
            u = pi*a*sin(theta)/lambda;
            I = (sin(u)./u).^2;
            axis(app.UIAxes,'xy');
            view(app.UIAxes,2);
            plot(app.UIAxes,x,I);
            % image(app.UIAxes,x,1,5000*I);              
            % colormap(app.UIAxes,gray);
            % axis(app.UIAxes,'on');
        end

        % Button pushed function: df_2d_button
        function df_2d_buttonPushed(app, event)
            global N;
            lambda = app.df_la_slider.Value*1e-9;
            a = app.df_a_slider.Value*1e-3;
            %a = (400:1:780)*1e-9;
            theta = linspace(-8e-5, 8e-5, N);  
            %theta = linspace(-5e-2, 5e-2, N);  
            x = tan(theta);
            u = pi*a*sin(theta)/lambda;
            I = (sin(u)./u).^2;
            view(app.UIAxes,2);
            colormap(app.UIAxes,gray);
            image(app.UIAxes,x,1,5000*I); 
            % axis(app.UIAxes,'off');
            
        end

        % Button pushed function: df_3d_button
        function df_3d_buttonButtonPushed(app, event)
            global N;
            lambda = app.df_la_slider.Value*1e-9;
            a = app.df_a_slider.Value*1e-3;
            theta = linspace(-8e-5, 8e-5, N);  
            %theta = linspace(-5e-2, 5e-2, N);  
            x = tan(theta);
            y = tan(theta);
            %y = linspace(-1, 1, 2);
            [X, Y] = meshgrid(x, y);
            alpha = pi*a*sin(theta)/lambda;
            I = (sin(alpha)./alpha).^2;     
            U=repmat(I,round(N),1);
            % axis(app.UIAxes,'xy');
            mesh(app.UIAxes,X, Y, U);
            colormap(app.UIAxes,'default');
        end

        % Button pushed function: df_dh_button
        function df_dh_buttonPushed(app, event)
            global N; global M;
            view(app.UIAxes,2);
            colormap(app.UIAxes,gray);
             switch app.df_DropDown.Value
                case "波长"
                    a = app.df_a_slider.Value*1e-3;
                    theta = linspace(-8e-5, 8e-5, N);  
                    x = tan(theta);      
                    for lambda = linspace(400,780,M)*1e-9
                        app.df_la_slider.Value = lambda*1e9;
                        app.df_la_spin.Value = lambda*1e9;
                        u = pi*a*sin(theta)/lambda;
                        I = (sin(u)./u).^2;
                        image(app.UIAxes,x,1,5000*I);  
                        %axis(app.UIAxes,'off');
                        drawnow;
                    end
                case "单缝宽度"
                    %a = 20e-3;
                    lambda = app.df_la_slider.Value*1e-9;
                    theta = linspace(-8e-5, 8e-5, N);  
                    x = tan(theta);      
                    for a = linspace(20,100,M)*1e-3
                        app.df_a_slider.Value = a*1e3;
                        app.df_a_spin.Value = a*1e3;
                        u = pi*a*sin(theta)/lambda;
                        I = (sin(u)./u).^2;            
                        image(app.UIAxes,x,1,5000*I);   
                        %axis(app.UIAxes,'off');
                        drawnow;
                    end
                 case "缝宽(多波长)"
                    lambda=[660,610,570,550,460,440,410]*1e-9; %七色
                    H = round(N); %图样高度
                    theta = linspace(-5e-5, 5e-5, N);
                    for a = linspace(20,100,M)*1e-3
                        app.df_a_slider.Value = a*1e3;
                        app.df_a_spin.Value = a*1e3;
                        RGB=zeros(length(lambda),3); %RGB结果矩阵置零
                        Iw=zeros(H,length(theta),3); %仿真光屏矩阵(仿真结果 RGB 
                        Irgb=zeros(H,length(theta),3); %用于记录各色光衍射结果的
                            for k=1:length(lambda)
                                alpha=pi*a*sin(theta)/lambda(k);
                                I=(sinc(alpha)).^2; %衍射的相对光强
                                RGB(k,:)=wave2rgb(lambda(k))./255;
                                Irgb(:,:,1)=repmat(I*RGB(k,1),H,1);
                                Irgb(:,:,2)=repmat(I*RGB(k,2),H,1);
                                Irgb(:,:,3)=repmat(I*RGB(k,3),H,1);
                                %计算白光光栅衍射 RGB 值图像矩阵数据
                                Iw=Iw+Irgb; %把各色光衍射的 RGB 值矩阵计入仿真结果 RGB 值图像矩阵中
                                Irgb=[];
                            end
                        Br=1/max(max(max(Iw))); %调整 Irgb 矩阵元素的最大值为 1 的系数
                        II=Iw*Br*100; %调节仿真图像亮度
                        % imshow(app.UIAxes,theta,theta,II); %显示仿真结果
                        imagesc(app.UIAxes,II);
                        % axis(app.UIAxes,'off');
                        drawnow;
                    end
             end
        end

        % Value changed function: duo_d_slider
        function duo_d_sliderValueChanged(app, event)
            app.duo_d_spin.Value = app.duo_d_slider.Value;
        end

        % Value changed function: duo_d_spin
        function duo_d_spinValueChanged(app, event)
            app.duo_d_slider.Value = app.duo_d_spin.Value;
        end

        % Value changed function: duo_n_slider
        function duo_n_sliderValueChanged(app, event)
            app.duo_n_spin.Value = app.duo_n_slider.Value;
        end

        % Value changed function: duo_n_spin
        function duo_n_spinValueChanged(app, event)
            app.duo_n_slider.Value = app.duo_n_spin.Value;
        end

        % Value changed function: duo_a_slider
        function duo_a_sliderValueChanged(app, event)
            app.duo_a_spin.Value = app.duo_a_slider.Value;
        end

        % Value changed function: duo_a_spin
        function duo_a_spinValueChanged(app, event)
            app.duo_a_slider.Value = app.duo_a_spin.Value;
        end

        % Value changed function: duo_la_slider
        function duo_la_sliderValueChanged(app, event)
            app.duo_la_spin.Value = app.duo_la_slider.Value;
        end

        % Value changed function: duo_la_spin
        function duo_la_spinValueChanged(app, event)
            app.duo_la_slider.Value = app.duo_la_spin.Value;
        end

        % Button pushed function: duo_ty_button
        function duo_ty_buttonPushed(app, event)
            global N;
            lambda = app.duo_la_slider.Value*1e-9;
            a = app.duo_a_slider.Value*1e-3;
            d = app.duo_d_slider.Value*1e-3;
            Ns = round(app.duo_n_slider.Value);
            %a = (400:1:780)*1e-9;
            theta = linspace(-8e-5, 8e-5, N);  
            %theta = linspace(-5e-2, 5e-2, N);  
            x = tan(theta);
            delta = 2*pi/lambda*d*sin(theta);
            alpha = pi*a*sin(theta)/lambda;
            I = (sin(alpha)./alpha).^2.*(sin(Ns.*delta./2)./sin(delta./2)).^2; 
            view(app.UIAxes,2);
            colormap(app.UIAxes,gray);
            image(app.UIAxes,x,1,5000*I);           
        end

        % Button pushed function: duo_1d_button
        function duo_1d_buttonButtonPushed(app, event)
            global N;
            lambda = app.duo_la_slider.Value*1e-9;
            a = app.duo_a_slider.Value*1e-3;
            d = app.duo_d_slider.Value*1e-3;
            Ns = round(app.duo_n_slider.Value);
            %a = (400:1:780)*1e-9;
            theta = linspace(-8e-5, 8e-5, N);  
            %theta = linspace(-5e-2, 5e-2, N);  
            x = tan(theta);
            delta = 2*pi/lambda*d*sin(theta);
            alpha = pi*a*sin(theta)/lambda;
            I = (sin(alpha)./alpha).^2.*(sin(Ns.*delta./2)./sin(delta./2)).^2; 
            view(app.UIAxes,2);
            axis(app.UIAxes,'xy');
            colormap(app.UIAxes,gray);
            %image(app.UIAxes,x,1,5000*I);   
            plot(app.UIAxes,x,I);
        end

        % Button pushed function: duo_dh_button
        function duo_2d_buttonPushed(app, event)
            global N; global M;
            view(app.UIAxes,2);
            colormap(app.UIAxes,gray);
            switch app.duo_DropDown.Value
                case "波长"
                    %lambda = app.duo_la_slider.Value*1e-9;
                    a = app.duo_a_slider.Value*1e-3;
                    d = app.duo_d_slider.Value*1e-3;
                    Ns = round(app.duo_n_slider.Value);
                    theta = linspace(-8e-5, 8e-5, N);  
                    x = tan(theta);
                    for lambda = linspace(400,780,M)*1e-9
                        app.duo_la_spin.Value = lambda*1e9;
                        app.duo_la_slider.Value = lambda*1e9;
                        delta = 2*pi/lambda*d*sin(theta);
                        alpha = pi*a*sin(theta)/lambda;
                        I = (sin(alpha)./alpha).^2.*(sin(Ns.*delta./2)./sin(delta./2)).^2; 
                        image(app.UIAxes,x,1,5000*I);   
                        drawnow;
                    end
                case "狭缝宽度"
                    lambda = app.duo_la_slider.Value*1e-9;
                    %a = app.duo_a_slider.Value*1e-3;
                    d = app.duo_d_slider.Value*1e-3;
                    Ns = round(app.duo_n_slider.Value);
                    theta = linspace(-8e-5, 8e-5, N);  
                    x = tan(theta);
                    for a = linspace(10,100,M)*1e-3
                        app.duo_a_spin.Value = a*1e3;
                        app.duo_a_slider.Value = a*1e3;
                        delta = 2*pi/lambda*d*sin(theta);
                        alpha = pi*a*sin(theta)/lambda;
                        I = (sin(alpha)./alpha).^2.*(sin(Ns.*delta./2)./sin(delta./2)).^2; 
                        image(app.UIAxes,x,1,5000*I);   
                        drawnow;
                    end
                case "缝间距"
                    lambda = app.duo_la_slider.Value*1e-9;
                    a = app.duo_a_slider.Value*1e-3;
                    %d = app.duo_d_slider.Value*1e-3;
                    Ns = round(app.duo_n_slider.Value);
                    theta = linspace(-8e-5, 8e-5, N);  
                    x = tan(theta);
                    for d = linspace(100,500,M)*1e-3
                        app.duo_d_spin.Value = d*1e3;
                        app.duo_d_slider.Value = d*1e3;
                        delta = 2*pi/lambda*d*sin(theta);
                        alpha = pi*a*sin(theta)/lambda;
                        I = (sin(alpha)./alpha).^2.*(sin(Ns.*delta./2)./sin(delta./2)).^2; 
                        image(app.UIAxes,x,1,5000*I);   
                        drawnow;
                    end
                case "狭缝数"
                    lambda = app.duo_la_slider.Value*1e-9;
                    a = app.duo_a_slider.Value*1e-3;
                    d = app.duo_d_slider.Value*1e-3;
                    %Ns = round(app.duo_n_slider.Value);
                    theta = linspace(-8e-5, 8e-5, N);  
                    x = tan(theta);
                    if M > 1
                        M = 11;
                    end
                    for Ns = round(linspace(1,11,M))
                        app.duo_n_spin.Value = Ns;
                        app.duo_n_slider.Value = Ns;
                        delta = 2*pi/lambda*d*sin(theta);
                        alpha = pi*a*sin(theta)/lambda;
                        I = (sin(alpha)./alpha).^2.*(sin(Ns.*delta./2)./sin(delta./2)).^2; 
                        image(app.UIAxes,x,1,5000*I);
                        drawnow;
                    end
            end
        end

        % Button pushed function: duo_3d_button
        function duo_3d_buttonButtonPushed(app, event)
            global N;
            lambda = app.duo_la_slider.Value*1e-9;
            a = app.duo_a_slider.Value*1e-3;
            d = app.duo_d_slider.Value*1e-3;
            Ns = round(app.duo_n_slider.Value);
            %a = (400:1:780)*1e-9;
            theta = linspace(-8e-5, 8e-5, N);  
            %theta = linspace(-5e-2, 5e-2, N);  
            x = tan(theta);
            y = tan(theta);
            [X, Y] = meshgrid(x, y);
            delta = 2*pi/lambda*d*sin(theta);
            alpha = pi*a*sin(theta)/lambda;
            I = (sin(alpha)./alpha).^2.*(sin(Ns.*delta./2)./sin(delta./2)).^2; 
            colormap(app.UIAxes,'default');
            %image(app.UIAxes,x,1,5000*I);   
            %plot(app.UIAxes,x,I);
            U=repmat(I,round(N),1);
            mesh(app.UIAxes,X,Y,U)
        end

        % Value changed function: jk_la_slider
        function jk_la_sliderValueChanged(app, event)
            app.jk_la_spin.Value = app.jk_la_slider.Value;
        end

        % Value changed function: jk_la_spin
        function jk_la_spinValueChanged(app, event)
            app.jk_la_slider.Value = app.jk_la_spin.Value;
        end

        % Value changed function: jk_a_slider
        function jk_a_sliderValueChanged(app, event)
            app.jk_a_spin.Value = app.jk_a_slider.Value;
        end

        % Value changed function: jk_a_spin
        function jk_a_spinValueChanged(app, event)
            app.jk_a_slider.Value = app.jk_a_spin.Value;
        end

        % Value changed function: jk_b_slider
        function jk_b_sliderValueChanged(app, event)
            app.jk_b_spin.Value = app.jk_b_slider.Value;
        end

        % Value changed function: jk_b_spin
        function jk_b_spinValueChanged(app, event)
            app.jk_b_slider.Value = app.jk_b_spin.Value;
        end

        % Button pushed function: jk_1d_button
        function jk_1d_buttonButtonPushed(app, event)
            global N;
            lambda = app.jk_la_slider.Value*1e-9;
            a = app.jk_a_slider.Value*1e-3;
            b = app.jk_b_slider.Value*1e-3;
            theta = linspace(-8e-5, 8e-5, N);  
            x = tan(theta);
            xx = 0.*theta;
            y = tan(theta);
            yy = 0.*theta;
            u1 = pi*a*sin(theta)/lambda;
            I1 = (sin(u1)./u1).^2;
            u2 = pi*b*sin(theta)/lambda;
            I2 = (sin(u2)./u2).^2;
            axis(app.UIAxes,'xy');
            view(app.UIAxes,2);
            plot3(app.UIAxes,x,xx,I1,yy,y,I2);
            % plot3(app.UIAxes,);
            %colormap(app.UIAxes,'default');
            %plot3(app.UIAxes,x,y,U(N/2,:))
            %imagesc(app.UIAxes,x, y, U.^0.35);
        end

        % Button pushed function: jk_2d_button
        function jk_2d_buttonButtonPushed(app, event)
            global N;
            lambda = app.jk_la_slider.Value*1e-9;
            a = app.jk_a_slider.Value*1e-3;
            b = app.jk_b_slider.Value*1e-3;
            theta = linspace(-8e-5, 8e-5, N);  
            x = tan(theta);
            y = tan(theta);
            alpha = pi*a*sin(theta)/lambda;
            I1 = (sin(alpha)./alpha).^2; 
            beta = pi*b*sin(theta)/lambda;
            I2 = (sin(beta)./beta).^2;    
            % N=200;
            U= zeros(N,N);
            for i=1:N
                for j=1:N
                    U(j,i)=I1(i)*I2(j);
                end
            end
            view(app.UIAxes,2);
            colormap(app.UIAxes,gray);
            imagesc(app.UIAxes,x, y, U.^0.35);
        end

        % Button pushed function: jk_3d_button
        function jk_3d_buttonButtonPushed(app, event)
            global N;
            lambda = app.jk_la_slider.Value*1e-9;
            a = app.jk_a_slider.Value*1e-3;
            b = app.jk_b_slider.Value*1e-3;
            theta = linspace(-8e-5, 8e-5, N);  
            x = tan(theta);
            y = tan(theta);
            [X, Y] = meshgrid(x, y);
            alpha = pi*a*sin(theta)/lambda;
            I1 = (sin(alpha)./alpha).^2; 
            beta = pi*b*sin(theta)/lambda;
            I2 = (sin(beta)./beta).^2;    
            % N=200;
            U= zeros(N,N);
            for i=1:N
                for j=1:N
                    U(j,i)=I1(i)*I2(j);
                end
            end
            colormap(app.UIAxes,'default');
            mesh(app.UIAxes,X,Y,U)
            %imagesc(app.UIAxes,x, y, U.^0.35);
        end

        % Button pushed function: jk_dh_button
        function jk_dh_buttonButtonPushed(app, event)
            global N; global M;
            view(app.UIAxes,2);
            colormap(app.UIAxes,gray);
            switch app.jk_DropDown.Value     
               case "波长"
                    %lambda = 632e-9;  % 光波长（632nm）
                    a = app.jk_a_slider.Value*1e-3;
                    b = app.jk_b_slider.Value*1e-3;
                    theta = linspace(-8e-5, 8e-5, N);  
                    x = tan(theta);
                    y = tan(theta);
                    for lambda = linspace(400,780,M)*1e-9
                        app.jk_la_spin.Value = lambda*1e9;
                        app.jk_la_slider.Value = lambda*1e9;
                        alpha = pi*a*sin(theta)/lambda;
                        I1 = (sin(alpha)./alpha).^2; 
                        beta = pi*b*sin(theta)/lambda;
                        I2 = (sin(beta)./beta).^2;  
                        % N=200;
                        U= zeros(N,N);
                        for i=1:N
                            for j=1:N
                                U(j,i)=I1(i)*I2(j);
                            end
                        end
                        imagesc(app.UIAxes,x, y, U.^0.35);    
                        
                       
                        drawnow;
                    end
             case "矩孔宽度"
                lambda = app.jk_la_slider.Value*1e-9;
                %a = app.jk_a_slider.Value*1e-3;
                b = app.jk_b_slider.Value*1e-3;
                theta = linspace(-8e-5, 8e-5, N);  
                x = tan(theta);
                y = tan(theta);
                for a = linspace(10,100,M)*1e-3
                    app.jk_a_spin.Value = a*1e3;
                    app.jk_a_slider.Value = a*1e3;
                    alpha = pi*a*sin(theta)/lambda;
                    I1 = (sin(alpha)./alpha).^2; 
                    beta = pi*b*sin(theta)/lambda;
                    I2 = (sin(beta)./beta).^2;     
                    % N=200;
                    U= zeros(N,N);
                    for i=1:N
                        for j=1:N
                            U(j,i)=I1(i)*I2(j);
                        end
                    end
                    imagesc(app.UIAxes,x, y, U.^0.35);
                    
                   
                    drawnow;
                end
              case "矩孔高度"
                lambda = app.jk_la_slider.Value*1e-9;
                a = app.jk_a_slider.Value*1e-3;
                %b = app.jk_b_slider.Value*1e-3;
                theta = linspace(-8e-5, 8e-5, N);  
                x = tan(theta);
                y = tan(theta);
                for b = linspace(10,100,M)*1e-3
                    app.jk_b_spin.Value = b*1e3;
                    app.jk_b_slider.Value = b*1e3;
                    alpha = pi*a*sin(theta)/lambda;
                    I1 = (sin(alpha)./alpha).^2; 
                    beta = pi*b*sin(theta)/lambda;
                    I2 = (sin(beta)./beta).^2;     
                    % N=200;
                    U= zeros(N,N);
                    for i=1:N
                        for j=1:N
                            U(j,i)=I1(i)*I2(j);
                        end
                    end
                    imagesc(app.UIAxes,x, y, U.^0.35);
                    drawnow;
                end
            end
        end

        % Value changed function: yk_la_slider
        function yk_la_sliderValueChanged(app, event)
            app.yk_la_spin.Value = app.yk_la_slider.Value;
        end

        % Value changed function: yk_la_spin
        function yk_la_spinValueChanged(app, event)
            app.yk_la_slider.Value = app.yk_la_spin.Value;
        end

        % Value changed function: yk_a_slider
        function yk_a_sliderValueChanged(app, event)
            app.yk_a_spin.Value = app.yk_a_slider.Value;
        end

        % Value changed function: yk_a_spin
        function yk_a_spinValueChanged(app, event)
            app.yk_a_slider.Value = app.yk_a_spin.Value;
        end

        % Button pushed function: yk_1d_button
        function yk_1d_buttonButtonPushed(app, event)
            global N;
            lambda = app.yk_la_slider.Value*1e-9;
            a = app.yk_a_slider.Value*1e-3;
            theta = linspace(-8e-5, 8e-5, N);  
            x = tan(theta);
            xx = 0.*theta;
            y = tan(theta);
            yy = 0.*theta;
            u = pi*a*sin(theta)/lambda;
            I = (sin(u)./u).^2;
            axis(app.UIAxes,'xy');
            view(app.UIAxes,2);
            plot3(app.UIAxes,x,xx,I,yy,y,I);
        end

        % Button pushed function: yk_2d_button
        function yk_2d_buttonButtonPushed(app, event)
            global N;
            lambda = app.yk_la_slider.Value*1e-9;
            a = app.yk_a_slider.Value*1e-3;
            %theta = linspace(-3e-5, 3e-5, N);  
            theta = linspace(-8e-5, 8e-5, N);  
            %f=50e-2;
            x = tan(theta);    
            y = tan(theta);    
            [X, Y] = meshgrid(x, y);
            Z = 2*pi*a*sqrt(X.^2+Y.^2)./(lambda);
            I = (2*besselj(1,Z)./Z).^2;
            view(app.UIAxes,2);
            colormap(app.UIAxes,gray);
            imagesc(app.UIAxes,x, y, I.^0.35);
            %image(app.UIAxes,x,1,5000*I);          
        end

        % Button pushed function: yk_3d_button
        function yk_3d_buttonButtonPushed(app, event)
            global N;
            lambda = app.yk_la_slider.Value*1e-9;
            a = app.yk_a_slider.Value*1e-3;
            %theta = linspace(-3e-5, 3e-5, N);  
            theta = linspace(-8e-5, 8e-5, N);  
            %f=50e-2;
            x = tan(theta);    
            y = tan(theta);               
            [X, Y] = meshgrid(x, y);
            Z = 2*pi*a*sqrt(X.^2+Y.^2)./(lambda);
            I = (2*besselj(1,Z)./Z).^2;
            colormap(app.UIAxes,'default');
            mesh(app.UIAxes,X,Y,I)

        end

        % Button pushed function: yk_dh_button
        function yk_dh_buttonButtonPushed(app, event)
            global N; global M;
            view(app.UIAxes,2);
            colormap(app.UIAxes,gray);
            switch app.yk_DropDown.Value
                case "波长"
                    %lambda = app.yk_la_slider.Value*1e-9;
                    a = app.yk_a_slider.Value*1e-3;
                    theta = linspace(-8e-5, 8e-5, N);  
                    x = tan(theta);    
                    y = tan(theta);    
                    [X, Y] = meshgrid(x, y);
                    for lambda = linspace(400,780,M)*1e-9
                        app.yk_la_spin.Value = lambda*1e9;
                        app.yk_la_slider.Value = lambda*1e9;
                        Z = 2*pi*a*sqrt(X.^2+Y.^2)./(lambda);
                        I = (2*besselj(1,Z)./Z).^2;
                        imagesc(app.UIAxes,x, y, I.^0.35);
                        %image(app.UIAxes,x,1,5000*I);
                        drawnow;
                    end
                case "圆孔半径"
                    lambda = app.yk_la_slider.Value*1e-9;
                    %a = app.yk_a_slider.Value*1e-3;
                    theta = linspace(-8e-5, 8e-5, N); 
                    x = tan(theta);    
                    y = tan(theta);    
                    [X, Y] = meshgrid(x, y);
                    for a = linspace(20,100,M)*1e-3
                        app.yk_a_spin.Value = a*1e3;
                        app.yk_a_slider.Value = a*1e3;
                        Z = 2*pi*a*sqrt(X.^2+Y.^2)./(lambda);
                        I = (2*besselj(1,Z)./Z).^2;
                        imagesc(app.UIAxes,x, y, I.^0.35);
                        %image(app.UIAxes,x,1,5000*I);      
                        drawnow;
                    end
            end
        end

        % Value changed function: set_p_slider
        function set_p_sliderValueChanged(app, event)
            app.set_p_spin.Value = app.set_p_slider.Value;
        end

        % Value changed function: set_p_spin
        function set_p_spinValueChanged(app, event)
             app.set_p_slider.Value = app.set_p_spin.Value;
        end

        % Value changed function: set_a_slider
        function set_a_sliderValueChanged(app, event)
            app.set_a_spin.Value = app.set_a_slider.Value;
        end

        % Value changed function: set_a_spin
        function set_a_spinValueChanged(app, event)
            app.set_a_slider.Value = app.set_a_spin.Value;
        end

        % Button pushed function: set_button
        function set_buttonButtonPushed(app, event)
            global N; global M;
            N = app.set_p_slider.Value;
            M = app.set_a_slider.Value;
        end

        % Button pushed function: set_defult_Button
        function set_defult_ButtonPushed(app, event)
            global N; global M;
            N = 500; M = 20;
            app.set_p_slider.Value = N;
            app.set_p_spin.Value = N;
            app.set_a_slider.Value = M;
            app.set_a_spin.Value = M;
        end

        % Value changed function: axis_switch
        function axis_switchValueChanged(app, event)
            if app.axis_switch.Value == '开'
                axis(app.UIAxes,'on');
            elseif app.axis_switch.Value == '关'
                axis(app.UIAxes,'off');
            end   
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.UIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {556, 556};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {345, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 894 556];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {345, '1x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create TabGroup
            app.TabGroup = uitabgroup(app.LeftPanel);
            app.TabGroup.Position = [1 6 338 549];

            % Create Tab
            app.Tab = uitab(app.TabGroup);
            app.Tab.Title = '单缝';

            % Create nmLabel
            app.nmLabel = uilabel(app.Tab);
            app.nmLabel.Position = [69 452 54 22];
            app.nmLabel.Text = '波长(nm)';

            % Create mLabel
            app.mLabel = uilabel(app.Tab);
            app.mLabel.Position = [190 452 79 22];
            app.mLabel.Text = '单缝宽度(μm)';

            % Create df_la_slider
            app.df_la_slider = uislider(app.Tab);
            app.df_la_slider.Limits = [400 780];
            app.df_la_slider.Orientation = 'vertical';
            app.df_la_slider.ValueChangedFcn = createCallbackFcn(app, @df_la_sliderValueChanged, true);
            app.df_la_slider.Position = [69 136 3 277];
            app.df_la_slider.Value = 632;

            % Create df_a_slider
            app.df_a_slider = uislider(app.Tab);
            app.df_a_slider.Limits = [10 100];
            app.df_a_slider.Orientation = 'vertical';
            app.df_a_slider.ValueChangedFcn = createCallbackFcn(app, @df_a_sliderValueChanged, true);
            app.df_a_slider.Position = [190 136 3 277];
            app.df_a_slider.Value = 50;

            % Create df_2d_button
            app.df_2d_button = uibutton(app.Tab, 'push');
            app.df_2d_button.ButtonPushedFcn = createCallbackFcn(app, @df_2d_buttonPushed, true);
            app.df_2d_button.Position = [120 94 75 26];
            app.df_2d_button.Text = '二维灰度图';

            % Create df_la_spin
            app.df_la_spin = uispinner(app.Tab);
            app.df_la_spin.Limits = [400 780];
            app.df_la_spin.RoundFractionalValues = 'on';
            app.df_la_spin.ValueDisplayFormat = '%d';
            app.df_la_spin.ValueChangedFcn = createCallbackFcn(app, @df_la_spinValueChanged, true);
            app.df_la_spin.HorizontalAlignment = 'center';
            app.df_la_spin.Position = [69 427 56 22];
            app.df_la_spin.Value = 632;

            % Create df_a_spin
            app.df_a_spin = uispinner(app.Tab);
            app.df_a_spin.Limits = [10 100];
            app.df_a_spin.RoundFractionalValues = 'on';
            app.df_a_spin.ValueDisplayFormat = '%d';
            app.df_a_spin.ValueChangedFcn = createCallbackFcn(app, @df_a_spinValueChanged, true);
            app.df_a_spin.HorizontalAlignment = 'center';
            app.df_a_spin.Position = [190 427 51 22];
            app.df_a_spin.Value = 50;

            % Create Label
            app.Label = uilabel(app.Tab);
            app.Label.HorizontalAlignment = 'right';
            app.Label.Position = [48 60 77 22];
            app.Label.Text = '变化的物理量';

            % Create df_DropDown
            app.df_DropDown = uidropdown(app.Tab);
            app.df_DropDown.Items = {'波长', '单缝宽度', '缝宽(多波长)'};
            app.df_DropDown.Position = [140 60 100 22];
            app.df_DropDown.Value = '波长';

            % Create df_dh_button
            app.df_dh_button = uibutton(app.Tab, 'push');
            app.df_dh_button.ButtonPushedFcn = createCallbackFcn(app, @df_dh_buttonPushed, true);
            app.df_dh_button.Position = [108 20 100 26];
            app.df_dh_button.Text = {'生成灰度图动画'; ''};

            % Create df_1d_button
            app.df_1d_button = uibutton(app.Tab, 'push');
            app.df_1d_button.ButtonPushedFcn = createCallbackFcn(app, @df_1d_buttonButtonPushed, true);
            app.df_1d_button.Position = [20 94 75 26];
            app.df_1d_button.Text = '强度曲线图';

            % Create df_3d_button
            app.df_3d_button = uibutton(app.Tab, 'push');
            app.df_3d_button.ButtonPushedFcn = createCallbackFcn(app, @df_3d_buttonButtonPushed, true);
            app.df_3d_button.Position = [220 94 75 26];
            app.df_3d_button.Text = '三维示意图';

            % Create Tab_2
            app.Tab_2 = uitab(app.TabGroup);
            app.Tab_2.Title = '多缝';

            % Create duo_dh_button
            app.duo_dh_button = uibutton(app.Tab_2, 'push');
            app.duo_dh_button.ButtonPushedFcn = createCallbackFcn(app, @duo_2d_buttonPushed, true);
            app.duo_dh_button.Position = [108 20 100 26];
            app.duo_dh_button.Text = {'生成灰度图动画'; ''};

            % Create nmLabel_2
            app.nmLabel_2 = uilabel(app.Tab_2);
            app.nmLabel_2.Position = [14 452 54 22];
            app.nmLabel_2.Text = '波长(nm)';

            % Create mLabel_2
            app.mLabel_2 = uilabel(app.Tab_2);
            app.mLabel_2.Position = [80 452 79 22];
            app.mLabel_2.Text = '狭缝宽度(μm)';

            % Create duo_la_slider
            app.duo_la_slider = uislider(app.Tab_2);
            app.duo_la_slider.Limits = [400 780];
            app.duo_la_slider.Orientation = 'vertical';
            app.duo_la_slider.ValueChangedFcn = createCallbackFcn(app, @duo_la_sliderValueChanged, true);
            app.duo_la_slider.Position = [14 136 3 277];
            app.duo_la_slider.Value = 632;

            % Create duo_a_slider
            app.duo_a_slider = uislider(app.Tab_2);
            app.duo_a_slider.Limits = [10 100];
            app.duo_a_slider.Orientation = 'vertical';
            app.duo_a_slider.ValueChangedFcn = createCallbackFcn(app, @duo_a_sliderValueChanged, true);
            app.duo_a_slider.Position = [91 136 3 277];
            app.duo_a_slider.Value = 50;

            % Create duo_ty_button
            app.duo_ty_button = uibutton(app.Tab_2, 'push');
            app.duo_ty_button.ButtonPushedFcn = createCallbackFcn(app, @duo_ty_buttonPushed, true);
            app.duo_ty_button.Position = [120 94 75 26];
            app.duo_ty_button.Text = {'二维灰度图'; ''};

            % Create duo_la_spin
            app.duo_la_spin = uispinner(app.Tab_2);
            app.duo_la_spin.Limits = [400 780];
            app.duo_la_spin.RoundFractionalValues = 'on';
            app.duo_la_spin.ValueDisplayFormat = '%d';
            app.duo_la_spin.ValueChangedFcn = createCallbackFcn(app, @duo_la_spinValueChanged, true);
            app.duo_la_spin.HorizontalAlignment = 'center';
            app.duo_la_spin.Position = [14 427 56 22];
            app.duo_la_spin.Value = 632;

            % Create duo_a_spin
            app.duo_a_spin = uispinner(app.Tab_2);
            app.duo_a_spin.Limits = [10 100];
            app.duo_a_spin.RoundFractionalValues = 'on';
            app.duo_a_spin.ValueDisplayFormat = '%d';
            app.duo_a_spin.ValueChangedFcn = createCallbackFcn(app, @duo_a_spinValueChanged, true);
            app.duo_a_spin.HorizontalAlignment = 'center';
            app.duo_a_spin.Position = [91 427 51 22];
            app.duo_a_spin.Value = 50;

            % Create DropDown_2Label
            app.DropDown_2Label = uilabel(app.Tab_2);
            app.DropDown_2Label.HorizontalAlignment = 'right';
            app.DropDown_2Label.Position = [48 60 77 22];
            app.DropDown_2Label.Text = '变化的物理量';

            % Create duo_DropDown
            app.duo_DropDown = uidropdown(app.Tab_2);
            app.duo_DropDown.Items = {'波长', '狭缝宽度', '缝间距', '狭缝数'};
            app.duo_DropDown.Position = [140 60 100 22];
            app.duo_DropDown.Value = '波长';

            % Create mmLabel
            app.mmLabel = uilabel(app.Tab_2);
            app.mmLabel.Position = [167 452 81 22];
            app.mmLabel.Text = '狭缝间距(mm)';

            % Create duo_d_slider
            app.duo_d_slider = uislider(app.Tab_2);
            app.duo_d_slider.Limits = [100 500];
            app.duo_d_slider.Orientation = 'vertical';
            app.duo_d_slider.ValueChangedFcn = createCallbackFcn(app, @duo_d_sliderValueChanged, true);
            app.duo_d_slider.Position = [174 136 3 277];
            app.duo_d_slider.Value = 200;

            % Create duo_d_spin
            app.duo_d_spin = uispinner(app.Tab_2);
            app.duo_d_spin.Limits = [100 500];
            app.duo_d_spin.RoundFractionalValues = 'on';
            app.duo_d_spin.ValueDisplayFormat = '%d';
            app.duo_d_spin.ValueChangedFcn = createCallbackFcn(app, @duo_d_spinValueChanged, true);
            app.duo_d_spin.HorizontalAlignment = 'center';
            app.duo_d_spin.Position = [173 427 58 22];
            app.duo_d_spin.Value = 200;

            % Create Label_2
            app.Label_2 = uilabel(app.Tab_2);
            app.Label_2.Position = [262 452 53 22];
            app.Label_2.Text = '狭缝数量';

            % Create duo_n_slider
            app.duo_n_slider = uislider(app.Tab_2);
            app.duo_n_slider.Limits = [1 11];
            app.duo_n_slider.Orientation = 'vertical';
            app.duo_n_slider.ValueChangedFcn = createCallbackFcn(app, @duo_n_sliderValueChanged, true);
            app.duo_n_slider.MinorTicks = [1 2 3 4 5 6 7 8 9 10 11];
            app.duo_n_slider.Position = [256 136 3 277];
            app.duo_n_slider.Value = 3;

            % Create duo_n_spin
            app.duo_n_spin = uispinner(app.Tab_2);
            app.duo_n_spin.Limits = [1 11];
            app.duo_n_spin.RoundFractionalValues = 'on';
            app.duo_n_spin.ValueDisplayFormat = '%d';
            app.duo_n_spin.ValueChangedFcn = createCallbackFcn(app, @duo_n_spinValueChanged, true);
            app.duo_n_spin.HorizontalAlignment = 'center';
            app.duo_n_spin.FontSize = 11;
            app.duo_n_spin.Position = [256 427 51 22];
            app.duo_n_spin.Value = 3;

            % Create duo_1d_button
            app.duo_1d_button = uibutton(app.Tab_2, 'push');
            app.duo_1d_button.ButtonPushedFcn = createCallbackFcn(app, @duo_1d_buttonButtonPushed, true);
            app.duo_1d_button.Position = [20 94 75 26];
            app.duo_1d_button.Text = '强度曲线图';

            % Create duo_3d_button
            app.duo_3d_button = uibutton(app.Tab_2, 'push');
            app.duo_3d_button.ButtonPushedFcn = createCallbackFcn(app, @duo_3d_buttonButtonPushed, true);
            app.duo_3d_button.Position = [220 94 75 26];
            app.duo_3d_button.Text = '三维示意图';

            % Create Tab_3
            app.Tab_3 = uitab(app.TabGroup);
            app.Tab_3.Title = '矩孔';

            % Create nmLabel_3
            app.nmLabel_3 = uilabel(app.Tab_3);
            app.nmLabel_3.Position = [29 452 54 22];
            app.nmLabel_3.Text = '波长(nm)';

            % Create mLabel_3
            app.mLabel_3 = uilabel(app.Tab_3);
            app.mLabel_3.Position = [115 452 79 22];
            app.mLabel_3.Text = '矩孔宽度(μm)';

            % Create jk_la_slider
            app.jk_la_slider = uislider(app.Tab_3);
            app.jk_la_slider.Limits = [400 780];
            app.jk_la_slider.Orientation = 'vertical';
            app.jk_la_slider.ValueChangedFcn = createCallbackFcn(app, @jk_la_sliderValueChanged, true);
            app.jk_la_slider.Position = [29 136 3 277];
            app.jk_la_slider.Value = 632;

            % Create jk_a_slider
            app.jk_a_slider = uislider(app.Tab_3);
            app.jk_a_slider.Limits = [10 100];
            app.jk_a_slider.Orientation = 'vertical';
            app.jk_a_slider.ValueChangedFcn = createCallbackFcn(app, @jk_a_sliderValueChanged, true);
            app.jk_a_slider.Position = [115 136 3 277];
            app.jk_a_slider.Value = 50;

            % Create jk_2d_button
            app.jk_2d_button = uibutton(app.Tab_3, 'push');
            app.jk_2d_button.ButtonPushedFcn = createCallbackFcn(app, @jk_2d_buttonButtonPushed, true);
            app.jk_2d_button.Position = [120 94 75 26];
            app.jk_2d_button.Text = '二维灰度图';

            % Create jk_la_spin
            app.jk_la_spin = uispinner(app.Tab_3);
            app.jk_la_spin.Limits = [400 780];
            app.jk_la_spin.RoundFractionalValues = 'on';
            app.jk_la_spin.ValueDisplayFormat = '%d';
            app.jk_la_spin.ValueChangedFcn = createCallbackFcn(app, @jk_la_spinValueChanged, true);
            app.jk_la_spin.HorizontalAlignment = 'center';
            app.jk_la_spin.Position = [29 427 56 22];
            app.jk_la_spin.Value = 632;

            % Create jk_a_spin
            app.jk_a_spin = uispinner(app.Tab_3);
            app.jk_a_spin.Limits = [10 100];
            app.jk_a_spin.RoundFractionalValues = 'on';
            app.jk_a_spin.ValueDisplayFormat = '%d';
            app.jk_a_spin.ValueChangedFcn = createCallbackFcn(app, @jk_a_spinValueChanged, true);
            app.jk_a_spin.HorizontalAlignment = 'center';
            app.jk_a_spin.Position = [115 427 51 22];
            app.jk_a_spin.Value = 50;

            % Create jk_dh_button
            app.jk_dh_button = uibutton(app.Tab_3, 'push');
            app.jk_dh_button.ButtonPushedFcn = createCallbackFcn(app, @jk_dh_buttonButtonPushed, true);
            app.jk_dh_button.Position = [108 20 100 26];
            app.jk_dh_button.Text = {'生成灰度图动画'; ''};

            % Create Label_3
            app.Label_3 = uilabel(app.Tab_3);
            app.Label_3.HorizontalAlignment = 'right';
            app.Label_3.Position = [48 60 77 22];
            app.Label_3.Text = '变化的物理量';

            % Create jk_DropDown
            app.jk_DropDown = uidropdown(app.Tab_3);
            app.jk_DropDown.Items = {'波长', '矩孔宽度', '矩孔高度'};
            app.jk_DropDown.Position = [140 60 100 22];
            app.jk_DropDown.Value = '波长';

            % Create mLabel_5
            app.mLabel_5 = uilabel(app.Tab_3);
            app.mLabel_5.Position = [207 452 79 22];
            app.mLabel_5.Text = '矩孔高度(μm)';

            % Create jk_b_slider
            app.jk_b_slider = uislider(app.Tab_3);
            app.jk_b_slider.Limits = [10 100];
            app.jk_b_slider.Orientation = 'vertical';
            app.jk_b_slider.ValueChangedFcn = createCallbackFcn(app, @jk_b_sliderValueChanged, true);
            app.jk_b_slider.Position = [207 136 3 277];
            app.jk_b_slider.Value = 50;

            % Create jk_b_spin
            app.jk_b_spin = uispinner(app.Tab_3);
            app.jk_b_spin.Limits = [10 100];
            app.jk_b_spin.RoundFractionalValues = 'on';
            app.jk_b_spin.ValueDisplayFormat = '%d';
            app.jk_b_spin.ValueChangedFcn = createCallbackFcn(app, @jk_b_spinValueChanged, true);
            app.jk_b_spin.HorizontalAlignment = 'center';
            app.jk_b_spin.Position = [207 427 51 22];
            app.jk_b_spin.Value = 50;

            % Create jk_1d_button
            app.jk_1d_button = uibutton(app.Tab_3, 'push');
            app.jk_1d_button.ButtonPushedFcn = createCallbackFcn(app, @jk_1d_buttonButtonPushed, true);
            app.jk_1d_button.Position = [20 94 75 26];
            app.jk_1d_button.Text = '强度曲线图';

            % Create jk_3d_button
            app.jk_3d_button = uibutton(app.Tab_3, 'push');
            app.jk_3d_button.ButtonPushedFcn = createCallbackFcn(app, @jk_3d_buttonButtonPushed, true);
            app.jk_3d_button.Position = [220 94 75 26];
            app.jk_3d_button.Text = '三维示意图';

            % Create Tab_4
            app.Tab_4 = uitab(app.TabGroup);
            app.Tab_4.Title = '圆孔';

            % Create nmLabel_4
            app.nmLabel_4 = uilabel(app.Tab_4);
            app.nmLabel_4.Position = [69 452 54 22];
            app.nmLabel_4.Text = '波长(nm)';

            % Create mLabel_4
            app.mLabel_4 = uilabel(app.Tab_4);
            app.mLabel_4.Position = [190 452 79 22];
            app.mLabel_4.Text = '圆孔半径(μm)';

            % Create yk_la_slider
            app.yk_la_slider = uislider(app.Tab_4);
            app.yk_la_slider.Limits = [400 780];
            app.yk_la_slider.Orientation = 'vertical';
            app.yk_la_slider.ValueChangedFcn = createCallbackFcn(app, @yk_la_sliderValueChanged, true);
            app.yk_la_slider.Position = [69 136 3 277];
            app.yk_la_slider.Value = 632;

            % Create yk_a_slider
            app.yk_a_slider = uislider(app.Tab_4);
            app.yk_a_slider.Limits = [10 100];
            app.yk_a_slider.Orientation = 'vertical';
            app.yk_a_slider.ValueChangedFcn = createCallbackFcn(app, @yk_a_sliderValueChanged, true);
            app.yk_a_slider.Position = [190 136 3 277];
            app.yk_a_slider.Value = 50;

            % Create yk_2d_button
            app.yk_2d_button = uibutton(app.Tab_4, 'push');
            app.yk_2d_button.ButtonPushedFcn = createCallbackFcn(app, @yk_2d_buttonButtonPushed, true);
            app.yk_2d_button.Position = [120 94 75 26];
            app.yk_2d_button.Text = '二维灰度图';

            % Create yk_la_spin
            app.yk_la_spin = uispinner(app.Tab_4);
            app.yk_la_spin.Limits = [400 780];
            app.yk_la_spin.RoundFractionalValues = 'on';
            app.yk_la_spin.ValueDisplayFormat = '%d';
            app.yk_la_spin.ValueChangedFcn = createCallbackFcn(app, @yk_la_spinValueChanged, true);
            app.yk_la_spin.HorizontalAlignment = 'center';
            app.yk_la_spin.Position = [69 427 56 22];
            app.yk_la_spin.Value = 632;

            % Create yk_a_spin
            app.yk_a_spin = uispinner(app.Tab_4);
            app.yk_a_spin.Limits = [10 100];
            app.yk_a_spin.RoundFractionalValues = 'on';
            app.yk_a_spin.ValueDisplayFormat = '%d';
            app.yk_a_spin.ValueChangedFcn = createCallbackFcn(app, @yk_a_spinValueChanged, true);
            app.yk_a_spin.HorizontalAlignment = 'center';
            app.yk_a_spin.Position = [190 427 51 22];
            app.yk_a_spin.Value = 50;

            % Create yk_dh_button
            app.yk_dh_button = uibutton(app.Tab_4, 'push');
            app.yk_dh_button.ButtonPushedFcn = createCallbackFcn(app, @yk_dh_buttonButtonPushed, true);
            app.yk_dh_button.Position = [108 20 100 26];
            app.yk_dh_button.Text = {'生成灰度图动画'; ''};

            % Create Label_4
            app.Label_4 = uilabel(app.Tab_4);
            app.Label_4.HorizontalAlignment = 'right';
            app.Label_4.Position = [48 60 77 22];
            app.Label_4.Text = '变化的物理量';

            % Create yk_DropDown
            app.yk_DropDown = uidropdown(app.Tab_4);
            app.yk_DropDown.Items = {'波长', '圆孔半径', ''};
            app.yk_DropDown.Position = [140 60 100 22];
            app.yk_DropDown.Value = '波长';

            % Create yk_1d_button
            app.yk_1d_button = uibutton(app.Tab_4, 'push');
            app.yk_1d_button.ButtonPushedFcn = createCallbackFcn(app, @yk_1d_buttonButtonPushed, true);
            app.yk_1d_button.Position = [20 94 75 26];
            app.yk_1d_button.Text = '强度曲线图';

            % Create yk_3d_button
            app.yk_3d_button = uibutton(app.Tab_4, 'push');
            app.yk_3d_button.ButtonPushedFcn = createCallbackFcn(app, @yk_3d_buttonButtonPushed, true);
            app.yk_3d_button.Position = [220 94 75 26];
            app.yk_3d_button.Text = '三维示意图';

            % Create Tab_5
            app.Tab_5 = uitab(app.TabGroup);
            app.Tab_5.Title = '设置';

            % Create Label_5
            app.Label_5 = uilabel(app.Tab_5);
            app.Label_5.Position = [42 452 65 22];
            app.Label_5.Text = '图样解析度';

            % Create Label_6
            app.Label_6 = uilabel(app.Tab_5);
            app.Label_6.Position = [142 452 53 22];
            app.Label_6.Text = '动画系数';

            % Create set_p_slider
            app.set_p_slider = uislider(app.Tab_5);
            app.set_p_slider.Limits = [100 2000];
            app.set_p_slider.Orientation = 'vertical';
            app.set_p_slider.ValueChangedFcn = createCallbackFcn(app, @set_p_sliderValueChanged, true);
            app.set_p_slider.Position = [42 136 3 277];
            app.set_p_slider.Value = 500;

            % Create set_a_slider
            app.set_a_slider = uislider(app.Tab_5);
            app.set_a_slider.Limits = [5 50];
            app.set_a_slider.Orientation = 'vertical';
            app.set_a_slider.ValueChangedFcn = createCallbackFcn(app, @set_a_sliderValueChanged, true);
            app.set_a_slider.Position = [142 136 3 277];
            app.set_a_slider.Value = 20;

            % Create set_p_spin
            app.set_p_spin = uispinner(app.Tab_5);
            app.set_p_spin.Limits = [100 2000];
            app.set_p_spin.RoundFractionalValues = 'on';
            app.set_p_spin.ValueDisplayFormat = '%d';
            app.set_p_spin.ValueChangedFcn = createCallbackFcn(app, @set_p_spinValueChanged, true);
            app.set_p_spin.HorizontalAlignment = 'center';
            app.set_p_spin.Position = [42 427 65 22];
            app.set_p_spin.Value = 500;

            % Create set_a_spin
            app.set_a_spin = uispinner(app.Tab_5);
            app.set_a_spin.Limits = [5 50];
            app.set_a_spin.RoundFractionalValues = 'on';
            app.set_a_spin.ValueDisplayFormat = '%d';
            app.set_a_spin.ValueChangedFcn = createCallbackFcn(app, @set_a_spinValueChanged, true);
            app.set_a_spin.HorizontalAlignment = 'center';
            app.set_a_spin.Position = [142 427 51 22];
            app.set_a_spin.Value = 20;

            % Create set_button
            app.set_button = uibutton(app.Tab_5, 'push');
            app.set_button.ButtonPushedFcn = createCallbackFcn(app, @set_buttonButtonPushed, true);
            app.set_button.Position = [108 11 100 26];
            app.set_button.Text = '保存';

            % Create set_Label
            app.set_Label = uilabel(app.Tab_5);
            app.set_Label.VerticalAlignment = 'top';
            app.set_Label.Position = [18 35 309 58];
            app.set_Label.Text = {'说明：'; '图样解析度越大，图样越清晰，过大可能导致运行缓慢；'; '动画系数越大，变量切换越慢，过大可能导致运行缓慢。'};

            % Create set_defult_Button
            app.set_defult_Button = uibutton(app.Tab_5, 'push');
            app.set_defult_Button.ButtonPushedFcn = createCallbackFcn(app, @set_defult_ButtonPushed, true);
            app.set_defult_Button.Position = [108 94 100 26];
            app.set_defult_Button.Text = '重置为默认值';

            % Create Label_7
            app.Label_7 = uilabel(app.Tab_5);
            app.Label_7.HorizontalAlignment = 'center';
            app.Label_7.Position = [233 258 65 22];
            app.Label_7.Text = '坐标轴显示';

            % Create axis_switch
            app.axis_switch = uiswitch(app.Tab_5, 'slider');
            app.axis_switch.Items = {'开', '关'};
            app.axis_switch.ValueChangedFcn = createCallbackFcn(app, @axis_switchValueChanged, true);
            app.axis_switch.Position = [242 283 45 20];
            app.axis_switch.Value = '开';

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;

            % Create UIAxes
            app.UIAxes = uiaxes(app.RightPanel);
            title(app.UIAxes, '夫琅禾费衍射仿真图样')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            app.UIAxes.TitleFontWeight = 'bold';
            app.UIAxes.Position = [1 6 542 549];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = myapp

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
endclassdef myapp < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure           matlab.ui.Figure
        GridLayout         matlab.ui.container.GridLayout
        LeftPanel          matlab.ui.container.Panel
        TabGroup           matlab.ui.container.TabGroup
        Tab                matlab.ui.container.Tab
        nmLabel            matlab.ui.control.Label
        mLabel             matlab.ui.control.Label
        df_la_slider       matlab.ui.control.Slider
        df_a_slider        matlab.ui.control.Slider
        df_2d_button       matlab.ui.control.Button
        df_la_spin         matlab.ui.control.Spinner
        df_a_spin          matlab.ui.control.Spinner
        Label              matlab.ui.control.Label
        df_DropDown        matlab.ui.control.DropDown
        df_dh_button       matlab.ui.control.Button
        df_1d_button       matlab.ui.control.Button
        df_3d_button       matlab.ui.control.Button
        Tab_2              matlab.ui.container.Tab
        duo_dh_button      matlab.ui.control.Button
        nmLabel_2          matlab.ui.control.Label
        mLabel_2           matlab.ui.control.Label
        duo_la_slider      matlab.ui.control.Slider
        duo_a_slider       matlab.ui.control.Slider
        duo_ty_button      matlab.ui.control.Button
        duo_la_spin        matlab.ui.control.Spinner
        duo_a_spin         matlab.ui.control.Spinner
        DropDown_2Label    matlab.ui.control.Label
        duo_DropDown       matlab.ui.control.DropDown
        mmLabel            matlab.ui.control.Label
        duo_d_slider       matlab.ui.control.Slider
        duo_d_spin         matlab.ui.control.Spinner
        Label_2            matlab.ui.control.Label
        duo_n_slider       matlab.ui.control.Slider
        duo_n_spin         matlab.ui.control.Spinner
        duo_1d_button      matlab.ui.control.Button
        duo_3d_button      matlab.ui.control.Button
        Tab_3              matlab.ui.container.Tab
        nmLabel_3          matlab.ui.control.Label
        mLabel_3           matlab.ui.control.Label
        jk_la_slider       matlab.ui.control.Slider
        jk_a_slider        matlab.ui.control.Slider
        jk_2d_button       matlab.ui.control.Button
        jk_la_spin         matlab.ui.control.Spinner
        jk_a_spin          matlab.ui.control.Spinner
        jk_dh_button       matlab.ui.control.Button
        Label_3            matlab.ui.control.Label
        jk_DropDown        matlab.ui.control.DropDown
        mLabel_5           matlab.ui.control.Label
        jk_b_slider        matlab.ui.control.Slider
        jk_b_spin          matlab.ui.control.Spinner
        jk_1d_button       matlab.ui.control.Button
        jk_3d_button       matlab.ui.control.Button
        Tab_4              matlab.ui.container.Tab
        nmLabel_4          matlab.ui.control.Label
        mLabel_4           matlab.ui.control.Label
        yk_la_slider       matlab.ui.control.Slider
        yk_a_slider        matlab.ui.control.Slider
        yk_2d_button       matlab.ui.control.Button
        yk_la_spin         matlab.ui.control.Spinner
        yk_a_spin          matlab.ui.control.Spinner
        yk_dh_button       matlab.ui.control.Button
        Label_4            matlab.ui.control.Label
        yk_DropDown        matlab.ui.control.DropDown
        yk_1d_button       matlab.ui.control.Button
        yk_3d_button       matlab.ui.control.Button
        Tab_5              matlab.ui.container.Tab
        Label_5            matlab.ui.control.Label
        Label_6            matlab.ui.control.Label
        set_p_slider       matlab.ui.control.Slider
        set_a_slider       matlab.ui.control.Slider
        set_p_spin         matlab.ui.control.Spinner
        set_a_spin         matlab.ui.control.Spinner
        set_button         matlab.ui.control.Button
        set_Label          matlab.ui.control.Label
        set_defult_Button  matlab.ui.control.Button
        Label_7            matlab.ui.control.Label
        axis_switch        matlab.ui.control.Switch
        RightPanel         matlab.ui.container.Panel
        UIAxes             matlab.ui.control.UIAxes
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.UIFigure.Name = "夫琅禾费衍射仿真工具";
            global N;
            global M;
            N = 500;
            M = 20;
            axis(app.UIAxes,'on');
            % axis(app.UIAxes,'off');
        end

        % Value changed function: df_la_slider
        function df_la_sliderValueChanged(app, event)
            app.df_la_spin.Value = app.df_la_slider.Value;
        end

        % Value changed function: df_la_spin
        function df_la_spinValueChanged(app, event)
            app.df_la_slider.Value = app.df_la_spin.Value;
        end

        % Value changed function: df_a_slider
        function df_a_sliderValueChanged(app, event)
            app.df_a_spin.Value = app.df_a_slider.Value;
        end

        % Value changed function: df_a_spin
        function df_a_spinValueChanged(app, event)
            app.df_a_slider.Value = app.df_a_spin.Value;
        end

        % Button pushed function: df_1d_button
        function df_1d_buttonButtonPushed(app, event)
            global N;
            lambda = app.df_la_slider.Value*1e-9;
            a = app.df_a_slider.Value*1e-3;
            %a = (400:1:780)*1e-9;
            theta = linspace(-8e-5, 8e-5, N);  
            %theta = linspace(-5e-2, 5e-2, N);  
            x = tan(theta);
            u = pi*a*sin(theta)/lambda;
            I = (sin(u)./u).^2;
            axis(app.UIAxes,'xy');
            view(app.UIAxes,2);
            plot(app.UIAxes,x,I);
            % image(app.UIAxes,x,1,5000*I);              
            % colormap(app.UIAxes,gray);
            % axis(app.UIAxes,'on');
        end

        % Button pushed function: df_2d_button
        function df_2d_buttonPushed(app, event)
            global N;
            lambda = app.df_la_slider.Value*1e-9;
            a = app.df_a_slider.Value*1e-3;
            %a = (400:1:780)*1e-9;
            theta = linspace(-8e-5, 8e-5, N);  
            %theta = linspace(-5e-2, 5e-2, N);  
            x = tan(theta);
            u = pi*a*sin(theta)/lambda;
            I = (sin(u)./u).^2;
            view(app.UIAxes,2);
            colormap(app.UIAxes,gray);
            image(app.UIAxes,x,1,5000*I); 
            % axis(app.UIAxes,'off');
            
        end

        % Button pushed function: df_3d_button
        function df_3d_buttonButtonPushed(app, event)
            global N;
            lambda = app.df_la_slider.Value*1e-9;
            a = app.df_a_slider.Value*1e-3;
            theta = linspace(-8e-5, 8e-5, N);  
            %theta = linspace(-5e-2, 5e-2, N);  
            x = tan(theta);
            y = tan(theta);
            %y = linspace(-1, 1, 2);
            [X, Y] = meshgrid(x, y);
            alpha = pi*a*sin(theta)/lambda;
            I = (sin(alpha)./alpha).^2;     
            U=repmat(I,round(N),1);
            % axis(app.UIAxes,'xy');
            mesh(app.UIAxes,X, Y, U);
            colormap(app.UIAxes,'default');
        end

        % Button pushed function: df_dh_button
        function df_dh_buttonPushed(app, event)
            global N; global M;
            view(app.UIAxes,2);
            colormap(app.UIAxes,gray);
             switch app.df_DropDown.Value
                case "波长"
                    a = app.df_a_slider.Value*1e-3;
                    theta = linspace(-8e-5, 8e-5, N);  
                    x = tan(theta);      
                    for lambda = linspace(400,780,M)*1e-9
                        app.df_la_slider.Value = lambda*1e9;
                        app.df_la_spin.Value = lambda*1e9;
                        u = pi*a*sin(theta)/lambda;
                        I = (sin(u)./u).^2;
                        image(app.UIAxes,x,1,5000*I);  
                        %axis(app.UIAxes,'off');
                        drawnow;
                    end
                case "单缝宽度"
                    %a = 20e-3;
                    lambda = app.df_la_slider.Value*1e-9;
                    theta = linspace(-8e-5, 8e-5, N);  
                    x = tan(theta);      
                    for a = linspace(20,100,M)*1e-3
                        app.df_a_slider.Value = a*1e3;
                        app.df_a_spin.Value = a*1e3;
                        u = pi*a*sin(theta)/lambda;
                        I = (sin(u)./u).^2;            
                        image(app.UIAxes,x,1,5000*I);   
                        %axis(app.UIAxes,'off');
                        drawnow;
                    end
                 case "缝宽(多波长)"
                    lambda=[660,610,570,550,460,440,410]*1e-9; %七色
                    H = round(N); %图样高度
                    theta = linspace(-5e-5, 5e-5, N);
                    for a = linspace(20,100,M)*1e-3
                        app.df_a_slider.Value = a*1e3;
                        app.df_a_spin.Value = a*1e3;
                        RGB=zeros(length(lambda),3); %RGB结果矩阵置零
                        Iw=zeros(H,length(theta),3); %仿真光屏矩阵(仿真结果 RGB 
                        Irgb=zeros(H,length(theta),3); %用于记录各色光衍射结果的
                            for k=1:length(lambda)
                                alpha=pi*a*sin(theta)/lambda(k);
                                I=(sinc(alpha)).^2; %衍射的相对光强
                                RGB(k,:)=wave2rgb(lambda(k))./255;
                                Irgb(:,:,1)=repmat(I*RGB(k,1),H,1);
                                Irgb(:,:,2)=repmat(I*RGB(k,2),H,1);
                                Irgb(:,:,3)=repmat(I*RGB(k,3),H,1);
                                %计算白光光栅衍射 RGB 值图像矩阵数据
                                Iw=Iw+Irgb; %把各色光衍射的 RGB 值矩阵计入仿真结果 RGB 值图像矩阵中
                                Irgb=[];
                            end
                        Br=1/max(max(max(Iw))); %调整 Irgb 矩阵元素的最大值为 1 的系数
                        II=Iw*Br*100; %调节仿真图像亮度
                        % imshow(app.UIAxes,theta,theta,II); %显示仿真结果
                        imagesc(app.UIAxes,II);
                        % axis(app.UIAxes,'off');
                        drawnow;
                    end
             end
        end

        % Value changed function: duo_d_slider
        function duo_d_sliderValueChanged(app, event)
            app.duo_d_spin.Value = app.duo_d_slider.Value;
        end

        % Value changed function: duo_d_spin
        function duo_d_spinValueChanged(app, event)
            app.duo_d_slider.Value = app.duo_d_spin.Value;
        end

        % Value changed function: duo_n_slider
        function duo_n_sliderValueChanged(app, event)
            app.duo_n_spin.Value = app.duo_n_slider.Value;
        end

        % Value changed function: duo_n_spin
        function duo_n_spinValueChanged(app, event)
            app.duo_n_slider.Value = app.duo_n_spin.Value;
        end

        % Value changed function: duo_a_slider
        function duo_a_sliderValueChanged(app, event)
            app.duo_a_spin.Value = app.duo_a_slider.Value;
        end

        % Value changed function: duo_a_spin
        function duo_a_spinValueChanged(app, event)
            app.duo_a_slider.Value = app.duo_a_spin.Value;
        end

        % Value changed function: duo_la_slider
        function duo_la_sliderValueChanged(app, event)
            app.duo_la_spin.Value = app.duo_la_slider.Value;
        end

        % Value changed function: duo_la_spin
        function duo_la_spinValueChanged(app, event)
            app.duo_la_slider.Value = app.duo_la_spin.Value;
        end

        % Button pushed function: duo_ty_button
        function duo_ty_buttonPushed(app, event)
            global N;
            lambda = app.duo_la_slider.Value*1e-9;
            a = app.duo_a_slider.Value*1e-3;
            d = app.duo_d_slider.Value*1e-3;
            Ns = round(app.duo_n_slider.Value);
            %a = (400:1:780)*1e-9;
            theta = linspace(-8e-5, 8e-5, N);  
            %theta = linspace(-5e-2, 5e-2, N);  
            x = tan(theta);
            delta = 2*pi/lambda*d*sin(theta);
            alpha = pi*a*sin(theta)/lambda;
            I = (sin(alpha)./alpha).^2.*(sin(Ns.*delta./2)./sin(delta./2)).^2; 
            view(app.UIAxes,2);
            colormap(app.UIAxes,gray);
            image(app.UIAxes,x,1,5000*I);           
        end

        % Button pushed function: duo_1d_button
        function duo_1d_buttonButtonPushed(app, event)
            global N;
            lambda = app.duo_la_slider.Value*1e-9;
            a = app.duo_a_slider.Value*1e-3;
            d = app.duo_d_slider.Value*1e-3;
            Ns = round(app.duo_n_slider.Value);
            %a = (400:1:780)*1e-9;
            theta = linspace(-8e-5, 8e-5, N);  
            %theta = linspace(-5e-2, 5e-2, N);  
            x = tan(theta);
            delta = 2*pi/lambda*d*sin(theta);
            alpha = pi*a*sin(theta)/lambda;
            I = (sin(alpha)./alpha).^2.*(sin(Ns.*delta./2)./sin(delta./2)).^2; 
            view(app.UIAxes,2);
            axis(app.UIAxes,'xy');
            colormap(app.UIAxes,gray);
            %image(app.UIAxes,x,1,5000*I);   
            plot(app.UIAxes,x,I);
        end

        % Button pushed function: duo_dh_button
        function duo_2d_buttonPushed(app, event)
            global N; global M;
            view(app.UIAxes,2);
            colormap(app.UIAxes,gray);
            switch app.duo_DropDown.Value
                case "波长"
                    %lambda = app.duo_la_slider.Value*1e-9;
                    a = app.duo_a_slider.Value*1e-3;
                    d = app.duo_d_slider.Value*1e-3;
                    Ns = round(app.duo_n_slider.Value);
                    theta = linspace(-8e-5, 8e-5, N);  
                    x = tan(theta);
                    for lambda = linspace(400,780,M)*1e-9
                        app.duo_la_spin.Value = lambda*1e9;
                        app.duo_la_slider.Value = lambda*1e9;
                        delta = 2*pi/lambda*d*sin(theta);
                        alpha = pi*a*sin(theta)/lambda;
                        I = (sin(alpha)./alpha).^2.*(sin(Ns.*delta./2)./sin(delta./2)).^2; 
                        image(app.UIAxes,x,1,5000*I);   
                        drawnow;
                    end
                case "狭缝宽度"
                    lambda = app.duo_la_slider.Value*1e-9;
                    %a = app.duo_a_slider.Value*1e-3;
                    d = app.duo_d_slider.Value*1e-3;
                    Ns = round(app.duo_n_slider.Value);
                    theta = linspace(-8e-5, 8e-5, N);  
                    x = tan(theta);
                    for a = linspace(10,100,M)*1e-3
                        app.duo_a_spin.Value = a*1e3;
                        app.duo_a_slider.Value = a*1e3;
                        delta = 2*pi/lambda*d*sin(theta);
                        alpha = pi*a*sin(theta)/lambda;
                        I = (sin(alpha)./alpha).^2.*(sin(Ns.*delta./2)./sin(delta./2)).^2; 
                        image(app.UIAxes,x,1,5000*I);   
                        drawnow;
                    end
                case "缝间距"
                    lambda = app.duo_la_slider.Value*1e-9;
                    a = app.duo_a_slider.Value*1e-3;
                    %d = app.duo_d_slider.Value*1e-3;
                    Ns = round(app.duo_n_slider.Value);
                    theta = linspace(-8e-5, 8e-5, N);  
                    x = tan(theta);
                    for d = linspace(100,500,M)*1e-3
                        app.duo_d_spin.Value = d*1e3;
                        app.duo_d_slider.Value = d*1e3;
                        delta = 2*pi/lambda*d*sin(theta);
                        alpha = pi*a*sin(theta)/lambda;
                        I = (sin(alpha)./alpha).^2.*(sin(Ns.*delta./2)./sin(delta./2)).^2; 
                        image(app.UIAxes,x,1,5000*I);   
                        drawnow;
                    end
                case "狭缝数"
                    lambda = app.duo_la_slider.Value*1e-9;
                    a = app.duo_a_slider.Value*1e-3;
                    d = app.duo_d_slider.Value*1e-3;
                    %Ns = round(app.duo_n_slider.Value);
                    theta = linspace(-8e-5, 8e-5, N);  
                    x = tan(theta);
                    if M > 1
                        M = 11;
                    end
                    for Ns = round(linspace(1,11,M))
                        app.duo_n_spin.Value = Ns;
                        app.duo_n_slider.Value = Ns;
                        delta = 2*pi/lambda*d*sin(theta);
                        alpha = pi*a*sin(theta)/lambda;
                        I = (sin(alpha)./alpha).^2.*(sin(Ns.*delta./2)./sin(delta./2)).^2; 
                        image(app.UIAxes,x,1,5000*I);
                        drawnow;
                    end
            end
        end

        % Button pushed function: duo_3d_button
        function duo_3d_buttonButtonPushed(app, event)
            global N;
            lambda = app.duo_la_slider.Value*1e-9;
            a = app.duo_a_slider.Value*1e-3;
            d = app.duo_d_slider.Value*1e-3;
            Ns = round(app.duo_n_slider.Value);
            %a = (400:1:780)*1e-9;
            theta = linspace(-8e-5, 8e-5, N);  
            %theta = linspace(-5e-2, 5e-2, N);  
            x = tan(theta);
            y = tan(theta);
            [X, Y] = meshgrid(x, y);
            delta = 2*pi/lambda*d*sin(theta);
            alpha = pi*a*sin(theta)/lambda;
            I = (sin(alpha)./alpha).^2.*(sin(Ns.*delta./2)./sin(delta./2)).^2; 
            colormap(app.UIAxes,'default');
            %image(app.UIAxes,x,1,5000*I);   
            %plot(app.UIAxes,x,I);
            U=repmat(I,round(N),1);
            mesh(app.UIAxes,X,Y,U)
        end

        % Value changed function: jk_la_slider
        function jk_la_sliderValueChanged(app, event)
            app.jk_la_spin.Value = app.jk_la_slider.Value;
        end

        % Value changed function: jk_la_spin
        function jk_la_spinValueChanged(app, event)
            app.jk_la_slider.Value = app.jk_la_spin.Value;
        end

        % Value changed function: jk_a_slider
        function jk_a_sliderValueChanged(app, event)
            app.jk_a_spin.Value = app.jk_a_slider.Value;
        end

        % Value changed function: jk_a_spin
        function jk_a_spinValueChanged(app, event)
            app.jk_a_slider.Value = app.jk_a_spin.Value;
        end

        % Value changed function: jk_b_slider
        function jk_b_sliderValueChanged(app, event)
            app.jk_b_spin.Value = app.jk_b_slider.Value;
        end

        % Value changed function: jk_b_spin
        function jk_b_spinValueChanged(app, event)
            app.jk_b_slider.Value = app.jk_b_spin.Value;
        end

        % Button pushed function: jk_1d_button
        function jk_1d_buttonButtonPushed(app, event)
            global N;
            lambda = app.jk_la_slider.Value*1e-9;
            a = app.jk_a_slider.Value*1e-3;
            b = app.jk_b_slider.Value*1e-3;
            theta = linspace(-8e-5, 8e-5, N);  
            x = tan(theta);
            xx = 0.*theta;
            y = tan(theta);
            yy = 0.*theta;
            u1 = pi*a*sin(theta)/lambda;
            I1 = (sin(u1)./u1).^2;
            u2 = pi*b*sin(theta)/lambda;
            I2 = (sin(u2)./u2).^2;
            axis(app.UIAxes,'xy');
            view(app.UIAxes,2);
            plot3(app.UIAxes,x,xx,I1,yy,y,I2);
            % plot3(app.UIAxes,);
            %colormap(app.UIAxes,'default');
            %plot3(app.UIAxes,x,y,U(N/2,:))
            %imagesc(app.UIAxes,x, y, U.^0.35);
        end

        % Button pushed function: jk_2d_button
        function jk_2d_buttonButtonPushed(app, event)
            global N;
            lambda = app.jk_la_slider.Value*1e-9;
            a = app.jk_a_slider.Value*1e-3;
            b = app.jk_b_slider.Value*1e-3;
            theta = linspace(-8e-5, 8e-5, N);  
            x = tan(theta);
            y = tan(theta);
            alpha = pi*a*sin(theta)/lambda;
            I1 = (sin(alpha)./alpha).^2; 
            beta = pi*b*sin(theta)/lambda;
            I2 = (sin(beta)./beta).^2;    
            % N=200;
            U= zeros(N,N);
            for i=1:N
                for j=1:N
                    U(j,i)=I1(i)*I2(j);
                end
            end
            view(app.UIAxes,2);
            colormap(app.UIAxes,gray);
            imagesc(app.UIAxes,x, y, U.^0.35);
        end

        % Button pushed function: jk_3d_button
        function jk_3d_buttonButtonPushed(app, event)
            global N;
            lambda = app.jk_la_slider.Value*1e-9;
            a = app.jk_a_slider.Value*1e-3;
            b = app.jk_b_slider.Value*1e-3;
            theta = linspace(-8e-5, 8e-5, N);  
            x = tan(theta);
            y = tan(theta);
            [X, Y] = meshgrid(x, y);
            alpha = pi*a*sin(theta)/lambda;
            I1 = (sin(alpha)./alpha).^2; 
            beta = pi*b*sin(theta)/lambda;
            I2 = (sin(beta)./beta).^2;    
            % N=200;
            U= zeros(N,N);
            for i=1:N
                for j=1:N
                    U(j,i)=I1(i)*I2(j);
                end
            end
            colormap(app.UIAxes,'default');
            mesh(app.UIAxes,X,Y,U)
            %imagesc(app.UIAxes,x, y, U.^0.35);
        end

        % Button pushed function: jk_dh_button
        function jk_dh_buttonButtonPushed(app, event)
            global N; global M;
            view(app.UIAxes,2);
            colormap(app.UIAxes,gray);
            switch app.jk_DropDown.Value     
               case "波长"
                    %lambda = 632e-9;  % 光波长（632nm）
                    a = app.jk_a_slider.Value*1e-3;
                    b = app.jk_b_slider.Value*1e-3;
                    theta = linspace(-8e-5, 8e-5, N);  
                    x = tan(theta);
                    y = tan(theta);
                    for lambda = linspace(400,780,M)*1e-9
                        app.jk_la_spin.Value = lambda*1e9;
                        app.jk_la_slider.Value = lambda*1e9;
                        alpha = pi*a*sin(theta)/lambda;
                        I1 = (sin(alpha)./alpha).^2; 
                        beta = pi*b*sin(theta)/lambda;
                        I2 = (sin(beta)./beta).^2;  
                        % N=200;
                        U= zeros(N,N);
                        for i=1:N
                            for j=1:N
                                U(j,i)=I1(i)*I2(j);
                            end
                        end
                        imagesc(app.UIAxes,x, y, U.^0.35);    
                        
                       
                        drawnow;
                    end
             case "矩孔宽度"
                lambda = app.jk_la_slider.Value*1e-9;
                %a = app.jk_a_slider.Value*1e-3;
                b = app.jk_b_slider.Value*1e-3;
                theta = linspace(-8e-5, 8e-5, N);  
                x = tan(theta);
                y = tan(theta);
                for a = linspace(10,100,M)*1e-3
                    app.jk_a_spin.Value = a*1e3;
                    app.jk_a_slider.Value = a*1e3;
                    alpha = pi*a*sin(theta)/lambda;
                    I1 = (sin(alpha)./alpha).^2; 
                    beta = pi*b*sin(theta)/lambda;
                    I2 = (sin(beta)./beta).^2;     
                    % N=200;
                    U= zeros(N,N);
                    for i=1:N
                        for j=1:N
                            U(j,i)=I1(i)*I2(j);
                        end
                    end
                    imagesc(app.UIAxes,x, y, U.^0.35);
                    
                   
                    drawnow;
                end
              case "矩孔高度"
                lambda = app.jk_la_slider.Value*1e-9;
                a = app.jk_a_slider.Value*1e-3;
                %b = app.jk_b_slider.Value*1e-3;
                theta = linspace(-8e-5, 8e-5, N);  
                x = tan(theta);
                y = tan(theta);
                for b = linspace(10,100,M)*1e-3
                    app.jk_b_spin.Value = b*1e3;
                    app.jk_b_slider.Value = b*1e3;
                    alpha = pi*a*sin(theta)/lambda;
                    I1 = (sin(alpha)./alpha).^2; 
                    beta = pi*b*sin(theta)/lambda;
                    I2 = (sin(beta)./beta).^2;     
                    % N=200;
                    U= zeros(N,N);
                    for i=1:N
                        for j=1:N
                            U(j,i)=I1(i)*I2(j);
                        end
                    end
                    imagesc(app.UIAxes,x, y, U.^0.35);
                    drawnow;
                end
            end
        end

        % Value changed function: yk_la_slider
        function yk_la_sliderValueChanged(app, event)
            app.yk_la_spin.Value = app.yk_la_slider.Value;
        end

        % Value changed function: yk_la_spin
        function yk_la_spinValueChanged(app, event)
            app.yk_la_slider.Value = app.yk_la_spin.Value;
        end

        % Value changed function: yk_a_slider
        function yk_a_sliderValueChanged(app, event)
            app.yk_a_spin.Value = app.yk_a_slider.Value;
        end

        % Value changed function: yk_a_spin
        function yk_a_spinValueChanged(app, event)
            app.yk_a_slider.Value = app.yk_a_spin.Value;
        end

        % Button pushed function: yk_1d_button
        function yk_1d_buttonButtonPushed(app, event)
            global N;
            lambda = app.yk_la_slider.Value*1e-9;
            a = app.yk_a_slider.Value*1e-3;
            theta = linspace(-8e-5, 8e-5, N);  
            x = tan(theta);
            xx = 0.*theta;
            y = tan(theta);
            yy = 0.*theta;
            u = pi*a*sin(theta)/lambda;
            I = (sin(u)./u).^2;
            axis(app.UIAxes,'xy');
            view(app.UIAxes,2);
            plot3(app.UIAxes,x,xx,I,yy,y,I);
        end

        % Button pushed function: yk_2d_button
        function yk_2d_buttonButtonPushed(app, event)
            global N;
            lambda = app.yk_la_slider.Value*1e-9;
            a = app.yk_a_slider.Value*1e-3;
            %theta = linspace(-3e-5, 3e-5, N);  
            theta = linspace(-8e-5, 8e-5, N);  
            %f=50e-2;
            x = tan(theta);    
            y = tan(theta);    
            [X, Y] = meshgrid(x, y);
            Z = 2*pi*a*sqrt(X.^2+Y.^2)./(lambda);
            I = (2*besselj(1,Z)./Z).^2;
            view(app.UIAxes,2);
            colormap(app.UIAxes,gray);
            imagesc(app.UIAxes,x, y, I.^0.35);
            %image(app.UIAxes,x,1,5000*I);          
        end

        % Button pushed function: yk_3d_button
        function yk_3d_buttonButtonPushed(app, event)
            global N;
            lambda = app.yk_la_slider.Value*1e-9;
            a = app.yk_a_slider.Value*1e-3;
            %theta = linspace(-3e-5, 3e-5, N);  
            theta = linspace(-8e-5, 8e-5, N);  
            %f=50e-2;
            x = tan(theta);    
            y = tan(theta);               
            [X, Y] = meshgrid(x, y);
            Z = 2*pi*a*sqrt(X.^2+Y.^2)./(lambda);
            I = (2*besselj(1,Z)./Z).^2;
            colormap(app.UIAxes,'default');
            mesh(app.UIAxes,X,Y,I)

        end

        % Button pushed function: yk_dh_button
        function yk_dh_buttonButtonPushed(app, event)
            global N; global M;
            view(app.UIAxes,2);
            colormap(app.UIAxes,gray);
            switch app.yk_DropDown.Value
                case "波长"
                    %lambda = app.yk_la_slider.Value*1e-9;
                    a = app.yk_a_slider.Value*1e-3;
                    theta = linspace(-8e-5, 8e-5, N);  
                    x = tan(theta);    
                    y = tan(theta);    
                    [X, Y] = meshgrid(x, y);
                    for lambda = linspace(400,780,M)*1e-9
                        app.yk_la_spin.Value = lambda*1e9;
                        app.yk_la_slider.Value = lambda*1e9;
                        Z = 2*pi*a*sqrt(X.^2+Y.^2)./(lambda);
                        I = (2*besselj(1,Z)./Z).^2;
                        imagesc(app.UIAxes,x, y, I.^0.35);
                        %image(app.UIAxes,x,1,5000*I);
                        drawnow;
                    end
                case "圆孔半径"
                    lambda = app.yk_la_slider.Value*1e-9;
                    %a = app.yk_a_slider.Value*1e-3;
                    theta = linspace(-8e-5, 8e-5, N); 
                    x = tan(theta);    
                    y = tan(theta);    
                    [X, Y] = meshgrid(x, y);
                    for a = linspace(20,100,M)*1e-3
                        app.yk_a_spin.Value = a*1e3;
                        app.yk_a_slider.Value = a*1e3;
                        Z = 2*pi*a*sqrt(X.^2+Y.^2)./(lambda);
                        I = (2*besselj(1,Z)./Z).^2;
                        imagesc(app.UIAxes,x, y, I.^0.35);
                        %image(app.UIAxes,x,1,5000*I);      
                        drawnow;
                    end
            end
        end

        % Value changed function: set_p_slider
        function set_p_sliderValueChanged(app, event)
            app.set_p_spin.Value = app.set_p_slider.Value;
        end

        % Value changed function: set_p_spin
        function set_p_spinValueChanged(app, event)
             app.set_p_slider.Value = app.set_p_spin.Value;
        end

        % Value changed function: set_a_slider
        function set_a_sliderValueChanged(app, event)
            app.set_a_spin.Value = app.set_a_slider.Value;
        end

        % Value changed function: set_a_spin
        function set_a_spinValueChanged(app, event)
            app.set_a_slider.Value = app.set_a_spin.Value;
        end

        % Button pushed function: set_button
        function set_buttonButtonPushed(app, event)
            global N; global M;
            N = app.set_p_slider.Value;
            M = app.set_a_slider.Value;
        end

        % Button pushed function: set_defult_Button
        function set_defult_ButtonPushed(app, event)
            global N; global M;
            N = 500; M = 20;
            app.set_p_slider.Value = N;
            app.set_p_spin.Value = N;
            app.set_a_slider.Value = M;
            app.set_a_spin.Value = M;
        end

        % Value changed function: axis_switch
        function axis_switchValueChanged(app, event)
            if app.axis_switch.Value == '开'
                axis(app.UIAxes,'on');
            elseif app.axis_switch.Value == '关'
                axis(app.UIAxes,'off');
            end   
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.UIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {556, 556};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {345, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 894 556];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {345, '1x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create TabGroup
            app.TabGroup = uitabgroup(app.LeftPanel);
            app.TabGroup.Position = [1 6 338 549];

            % Create Tab
            app.Tab = uitab(app.TabGroup);
            app.Tab.Title = '单缝';

            % Create nmLabel
            app.nmLabel = uilabel(app.Tab);
            app.nmLabel.Position = [69 452 54 22];
            app.nmLabel.Text = '波长(nm)';

            % Create mLabel
            app.mLabel = uilabel(app.Tab);
            app.mLabel.Position = [190 452 79 22];
            app.mLabel.Text = '单缝宽度(μm)';

            % Create df_la_slider
            app.df_la_slider = uislider(app.Tab);
            app.df_la_slider.Limits = [400 780];
            app.df_la_slider.Orientation = 'vertical';
            app.df_la_slider.ValueChangedFcn = createCallbackFcn(app, @df_la_sliderValueChanged, true);
            app.df_la_slider.Position = [69 136 3 277];
            app.df_la_slider.Value = 632;

            % Create df_a_slider
            app.df_a_slider = uislider(app.Tab);
            app.df_a_slider.Limits = [10 100];
            app.df_a_slider.Orientation = 'vertical';
            app.df_a_slider.ValueChangedFcn = createCallbackFcn(app, @df_a_sliderValueChanged, true);
            app.df_a_slider.Position = [190 136 3 277];
            app.df_a_slider.Value = 50;

            % Create df_2d_button
            app.df_2d_button = uibutton(app.Tab, 'push');
            app.df_2d_button.ButtonPushedFcn = createCallbackFcn(app, @df_2d_buttonPushed, true);
            app.df_2d_button.Position = [120 94 75 26];
            app.df_2d_button.Text = '二维灰度图';

            % Create df_la_spin
            app.df_la_spin = uispinner(app.Tab);
            app.df_la_spin.Limits = [400 780];
            app.df_la_spin.RoundFractionalValues = 'on';
            app.df_la_spin.ValueDisplayFormat = '%d';
            app.df_la_spin.ValueChangedFcn = createCallbackFcn(app, @df_la_spinValueChanged, true);
            app.df_la_spin.HorizontalAlignment = 'center';
            app.df_la_spin.Position = [69 427 56 22];
            app.df_la_spin.Value = 632;

            % Create df_a_spin
            app.df_a_spin = uispinner(app.Tab);
            app.df_a_spin.Limits = [10 100];
            app.df_a_spin.RoundFractionalValues = 'on';
            app.df_a_spin.ValueDisplayFormat = '%d';
            app.df_a_spin.ValueChangedFcn = createCallbackFcn(app, @df_a_spinValueChanged, true);
            app.df_a_spin.HorizontalAlignment = 'center';
            app.df_a_spin.Position = [190 427 51 22];
            app.df_a_spin.Value = 50;

            % Create Label
            app.Label = uilabel(app.Tab);
            app.Label.HorizontalAlignment = 'right';
            app.Label.Position = [48 60 77 22];
            app.Label.Text = '变化的物理量';

            % Create df_DropDown
            app.df_DropDown = uidropdown(app.Tab);
            app.df_DropDown.Items = {'波长', '单缝宽度', '缝宽(多波长)'};
            app.df_DropDown.Position = [140 60 100 22];
            app.df_DropDown.Value = '波长';

            % Create df_dh_button
            app.df_dh_button = uibutton(app.Tab, 'push');
            app.df_dh_button.ButtonPushedFcn = createCallbackFcn(app, @df_dh_buttonPushed, true);
            app.df_dh_button.Position = [108 20 100 26];
            app.df_dh_button.Text = {'生成灰度图动画'; ''};

            % Create df_1d_button
            app.df_1d_button = uibutton(app.Tab, 'push');
            app.df_1d_button.ButtonPushedFcn = createCallbackFcn(app, @df_1d_buttonButtonPushed, true);
            app.df_1d_button.Position = [20 94 75 26];
            app.df_1d_button.Text = '强度曲线图';

            % Create df_3d_button
            app.df_3d_button = uibutton(app.Tab, 'push');
            app.df_3d_button.ButtonPushedFcn = createCallbackFcn(app, @df_3d_buttonButtonPushed, true);
            app.df_3d_button.Position = [220 94 75 26];
            app.df_3d_button.Text = '三维示意图';

            % Create Tab_2
            app.Tab_2 = uitab(app.TabGroup);
            app.Tab_2.Title = '多缝';

            % Create duo_dh_button
            app.duo_dh_button = uibutton(app.Tab_2, 'push');
            app.duo_dh_button.ButtonPushedFcn = createCallbackFcn(app, @duo_2d_buttonPushed, true);
            app.duo_dh_button.Position = [108 20 100 26];
            app.duo_dh_button.Text = {'生成灰度图动画'; ''};

            % Create nmLabel_2
            app.nmLabel_2 = uilabel(app.Tab_2);
            app.nmLabel_2.Position = [14 452 54 22];
            app.nmLabel_2.Text = '波长(nm)';

            % Create mLabel_2
            app.mLabel_2 = uilabel(app.Tab_2);
            app.mLabel_2.Position = [80 452 79 22];
            app.mLabel_2.Text = '狭缝宽度(μm)';

            % Create duo_la_slider
            app.duo_la_slider = uislider(app.Tab_2);
            app.duo_la_slider.Limits = [400 780];
            app.duo_la_slider.Orientation = 'vertical';
            app.duo_la_slider.ValueChangedFcn = createCallbackFcn(app, @duo_la_sliderValueChanged, true);
            app.duo_la_slider.Position = [14 136 3 277];
            app.duo_la_slider.Value = 632;

            % Create duo_a_slider
            app.duo_a_slider = uislider(app.Tab_2);
            app.duo_a_slider.Limits = [10 100];
            app.duo_a_slider.Orientation = 'vertical';
            app.duo_a_slider.ValueChangedFcn = createCallbackFcn(app, @duo_a_sliderValueChanged, true);
            app.duo_a_slider.Position = [91 136 3 277];
            app.duo_a_slider.Value = 50;

            % Create duo_ty_button
            app.duo_ty_button = uibutton(app.Tab_2, 'push');
            app.duo_ty_button.ButtonPushedFcn = createCallbackFcn(app, @duo_ty_buttonPushed, true);
            app.duo_ty_button.Position = [120 94 75 26];
            app.duo_ty_button.Text = {'二维灰度图'; ''};

            % Create duo_la_spin
            app.duo_la_spin = uispinner(app.Tab_2);
            app.duo_la_spin.Limits = [400 780];
            app.duo_la_spin.RoundFractionalValues = 'on';
            app.duo_la_spin.ValueDisplayFormat = '%d';
            app.duo_la_spin.ValueChangedFcn = createCallbackFcn(app, @duo_la_spinValueChanged, true);
            app.duo_la_spin.HorizontalAlignment = 'center';
            app.duo_la_spin.Position = [14 427 56 22];
            app.duo_la_spin.Value = 632;

            % Create duo_a_spin
            app.duo_a_spin = uispinner(app.Tab_2);
            app.duo_a_spin.Limits = [10 100];
            app.duo_a_spin.RoundFractionalValues = 'on';
            app.duo_a_spin.ValueDisplayFormat = '%d';
            app.duo_a_spin.ValueChangedFcn = createCallbackFcn(app, @duo_a_spinValueChanged, true);
            app.duo_a_spin.HorizontalAlignment = 'center';
            app.duo_a_spin.Position = [91 427 51 22];
            app.duo_a_spin.Value = 50;

            % Create DropDown_2Label
            app.DropDown_2Label = uilabel(app.Tab_2);
            app.DropDown_2Label.HorizontalAlignment = 'right';
            app.DropDown_2Label.Position = [48 60 77 22];
            app.DropDown_2Label.Text = '变化的物理量';

            % Create duo_DropDown
            app.duo_DropDown = uidropdown(app.Tab_2);
            app.duo_DropDown.Items = {'波长', '狭缝宽度', '缝间距', '狭缝数'};
            app.duo_DropDown.Position = [140 60 100 22];
            app.duo_DropDown.Value = '波长';

            % Create mmLabel
            app.mmLabel = uilabel(app.Tab_2);
            app.mmLabel.Position = [167 452 81 22];
            app.mmLabel.Text = '狭缝间距(mm)';

            % Create duo_d_slider
            app.duo_d_slider = uislider(app.Tab_2);
            app.duo_d_slider.Limits = [100 500];
            app.duo_d_slider.Orientation = 'vertical';
            app.duo_d_slider.ValueChangedFcn = createCallbackFcn(app, @duo_d_sliderValueChanged, true);
            app.duo_d_slider.Position = [174 136 3 277];
            app.duo_d_slider.Value = 200;

            % Create duo_d_spin
            app.duo_d_spin = uispinner(app.Tab_2);
            app.duo_d_spin.Limits = [100 500];
            app.duo_d_spin.RoundFractionalValues = 'on';
            app.duo_d_spin.ValueDisplayFormat = '%d';
            app.duo_d_spin.ValueChangedFcn = createCallbackFcn(app, @duo_d_spinValueChanged, true);
            app.duo_d_spin.HorizontalAlignment = 'center';
            app.duo_d_spin.Position = [173 427 58 22];
            app.duo_d_spin.Value = 200;

            % Create Label_2
            app.Label_2 = uilabel(app.Tab_2);
            app.Label_2.Position = [262 452 53 22];
            app.Label_2.Text = '狭缝数量';

            % Create duo_n_slider
            app.duo_n_slider = uislider(app.Tab_2);
            app.duo_n_slider.Limits = [1 11];
            app.duo_n_slider.Orientation = 'vertical';
            app.duo_n_slider.ValueChangedFcn = createCallbackFcn(app, @duo_n_sliderValueChanged, true);
            app.duo_n_slider.MinorTicks = [1 2 3 4 5 6 7 8 9 10 11];
            app.duo_n_slider.Position = [256 136 3 277];
            app.duo_n_slider.Value = 3;

            % Create duo_n_spin
            app.duo_n_spin = uispinner(app.Tab_2);
            app.duo_n_spin.Limits = [1 11];
            app.duo_n_spin.RoundFractionalValues = 'on';
            app.duo_n_spin.ValueDisplayFormat = '%d';
            app.duo_n_spin.ValueChangedFcn = createCallbackFcn(app, @duo_n_spinValueChanged, true);
            app.duo_n_spin.HorizontalAlignment = 'center';
            app.duo_n_spin.FontSize = 11;
            app.duo_n_spin.Position = [256 427 51 22];
            app.duo_n_spin.Value = 3;

            % Create duo_1d_button
            app.duo_1d_button = uibutton(app.Tab_2, 'push');
            app.duo_1d_button.ButtonPushedFcn = createCallbackFcn(app, @duo_1d_buttonButtonPushed, true);
            app.duo_1d_button.Position = [20 94 75 26];
            app.duo_1d_button.Text = '强度曲线图';

            % Create duo_3d_button
            app.duo_3d_button = uibutton(app.Tab_2, 'push');
            app.duo_3d_button.ButtonPushedFcn = createCallbackFcn(app, @duo_3d_buttonButtonPushed, true);
            app.duo_3d_button.Position = [220 94 75 26];
            app.duo_3d_button.Text = '三维示意图';

            % Create Tab_3
            app.Tab_3 = uitab(app.TabGroup);
            app.Tab_3.Title = '矩孔';

            % Create nmLabel_3
            app.nmLabel_3 = uilabel(app.Tab_3);
            app.nmLabel_3.Position = [29 452 54 22];
            app.nmLabel_3.Text = '波长(nm)';

            % Create mLabel_3
            app.mLabel_3 = uilabel(app.Tab_3);
            app.mLabel_3.Position = [115 452 79 22];
            app.mLabel_3.Text = '矩孔宽度(μm)';

            % Create jk_la_slider
            app.jk_la_slider = uislider(app.Tab_3);
            app.jk_la_slider.Limits = [400 780];
            app.jk_la_slider.Orientation = 'vertical';
            app.jk_la_slider.ValueChangedFcn = createCallbackFcn(app, @jk_la_sliderValueChanged, true);
            app.jk_la_slider.Position = [29 136 3 277];
            app.jk_la_slider.Value = 632;

            % Create jk_a_slider
            app.jk_a_slider = uislider(app.Tab_3);
            app.jk_a_slider.Limits = [10 100];
            app.jk_a_slider.Orientation = 'vertical';
            app.jk_a_slider.ValueChangedFcn = createCallbackFcn(app, @jk_a_sliderValueChanged, true);
            app.jk_a_slider.Position = [115 136 3 277];
            app.jk_a_slider.Value = 50;

            % Create jk_2d_button
            app.jk_2d_button = uibutton(app.Tab_3, 'push');
            app.jk_2d_button.ButtonPushedFcn = createCallbackFcn(app, @jk_2d_buttonButtonPushed, true);
            app.jk_2d_button.Position = [120 94 75 26];
            app.jk_2d_button.Text = '二维灰度图';

            % Create jk_la_spin
            app.jk_la_spin = uispinner(app.Tab_3);
            app.jk_la_spin.Limits = [400 780];
            app.jk_la_spin.RoundFractionalValues = 'on';
            app.jk_la_spin.ValueDisplayFormat = '%d';
            app.jk_la_spin.ValueChangedFcn = createCallbackFcn(app, @jk_la_spinValueChanged, true);
            app.jk_la_spin.HorizontalAlignment = 'center';
            app.jk_la_spin.Position = [29 427 56 22];
            app.jk_la_spin.Value = 632;

            % Create jk_a_spin
            app.jk_a_spin = uispinner(app.Tab_3);
            app.jk_a_spin.Limits = [10 100];
            app.jk_a_spin.RoundFractionalValues = 'on';
            app.jk_a_spin.ValueDisplayFormat = '%d';
            app.jk_a_spin.ValueChangedFcn = createCallbackFcn(app, @jk_a_spinValueChanged, true);
            app.jk_a_spin.HorizontalAlignment = 'center';
            app.jk_a_spin.Position = [115 427 51 22];
            app.jk_a_spin.Value = 50;

            % Create jk_dh_button
            app.jk_dh_button = uibutton(app.Tab_3, 'push');
            app.jk_dh_button.ButtonPushedFcn = createCallbackFcn(app, @jk_dh_buttonButtonPushed, true);
            app.jk_dh_button.Position = [108 20 100 26];
            app.jk_dh_button.Text = {'生成灰度图动画'; ''};

            % Create Label_3
            app.Label_3 = uilabel(app.Tab_3);
            app.Label_3.HorizontalAlignment = 'right';
            app.Label_3.Position = [48 60 77 22];
            app.Label_3.Text = '变化的物理量';

            % Create jk_DropDown
            app.jk_DropDown = uidropdown(app.Tab_3);
            app.jk_DropDown.Items = {'波长', '矩孔宽度', '矩孔高度'};
            app.jk_DropDown.Position = [140 60 100 22];
            app.jk_DropDown.Value = '波长';

            % Create mLabel_5
            app.mLabel_5 = uilabel(app.Tab_3);
            app.mLabel_5.Position = [207 452 79 22];
            app.mLabel_5.Text = '矩孔高度(μm)';

            % Create jk_b_slider
            app.jk_b_slider = uislider(app.Tab_3);
            app.jk_b_slider.Limits = [10 100];
            app.jk_b_slider.Orientation = 'vertical';
            app.jk_b_slider.ValueChangedFcn = createCallbackFcn(app, @jk_b_sliderValueChanged, true);
            app.jk_b_slider.Position = [207 136 3 277];
            app.jk_b_slider.Value = 50;

            % Create jk_b_spin
            app.jk_b_spin = uispinner(app.Tab_3);
            app.jk_b_spin.Limits = [10 100];
            app.jk_b_spin.RoundFractionalValues = 'on';
            app.jk_b_spin.ValueDisplayFormat = '%d';
            app.jk_b_spin.ValueChangedFcn = createCallbackFcn(app, @jk_b_spinValueChanged, true);
            app.jk_b_spin.HorizontalAlignment = 'center';
            app.jk_b_spin.Position = [207 427 51 22];
            app.jk_b_spin.Value = 50;

            % Create jk_1d_button
            app.jk_1d_button = uibutton(app.Tab_3, 'push');
            app.jk_1d_button.ButtonPushedFcn = createCallbackFcn(app, @jk_1d_buttonButtonPushed, true);
            app.jk_1d_button.Position = [20 94 75 26];
            app.jk_1d_button.Text = '强度曲线图';

            % Create jk_3d_button
            app.jk_3d_button = uibutton(app.Tab_3, 'push');
            app.jk_3d_button.ButtonPushedFcn = createCallbackFcn(app, @jk_3d_buttonButtonPushed, true);
            app.jk_3d_button.Position = [220 94 75 26];
            app.jk_3d_button.Text = '三维示意图';

            % Create Tab_4
            app.Tab_4 = uitab(app.TabGroup);
            app.Tab_4.Title = '圆孔';

            % Create nmLabel_4
            app.nmLabel_4 = uilabel(app.Tab_4);
            app.nmLabel_4.Position = [69 452 54 22];
            app.nmLabel_4.Text = '波长(nm)';

            % Create mLabel_4
            app.mLabel_4 = uilabel(app.Tab_4);
            app.mLabel_4.Position = [190 452 79 22];
            app.mLabel_4.Text = '圆孔半径(μm)';

            % Create yk_la_slider
            app.yk_la_slider = uislider(app.Tab_4);
            app.yk_la_slider.Limits = [400 780];
            app.yk_la_slider.Orientation = 'vertical';
            app.yk_la_slider.ValueChangedFcn = createCallbackFcn(app, @yk_la_sliderValueChanged, true);
            app.yk_la_slider.Position = [69 136 3 277];
            app.yk_la_slider.Value = 632;

            % Create yk_a_slider
            app.yk_a_slider = uislider(app.Tab_4);
            app.yk_a_slider.Limits = [10 100];
            app.yk_a_slider.Orientation = 'vertical';
            app.yk_a_slider.ValueChangedFcn = createCallbackFcn(app, @yk_a_sliderValueChanged, true);
            app.yk_a_slider.Position = [190 136 3 277];
            app.yk_a_slider.Value = 50;

            % Create yk_2d_button
            app.yk_2d_button = uibutton(app.Tab_4, 'push');
            app.yk_2d_button.ButtonPushedFcn = createCallbackFcn(app, @yk_2d_buttonButtonPushed, true);
            app.yk_2d_button.Position = [120 94 75 26];
            app.yk_2d_button.Text = '二维灰度图';

            % Create yk_la_spin
            app.yk_la_spin = uispinner(app.Tab_4);
            app.yk_la_spin.Limits = [400 780];
            app.yk_la_spin.RoundFractionalValues = 'on';
            app.yk_la_spin.ValueDisplayFormat = '%d';
            app.yk_la_spin.ValueChangedFcn = createCallbackFcn(app, @yk_la_spinValueChanged, true);
            app.yk_la_spin.HorizontalAlignment = 'center';
            app.yk_la_spin.Position = [69 427 56 22];
            app.yk_la_spin.Value = 632;

            % Create yk_a_spin
            app.yk_a_spin = uispinner(app.Tab_4);
            app.yk_a_spin.Limits = [10 100];
            app.yk_a_spin.RoundFractionalValues = 'on';
            app.yk_a_spin.ValueDisplayFormat = '%d';
            app.yk_a_spin.ValueChangedFcn = createCallbackFcn(app, @yk_a_spinValueChanged, true);
            app.yk_a_spin.HorizontalAlignment = 'center';
            app.yk_a_spin.Position = [190 427 51 22];
            app.yk_a_spin.Value = 50;

            % Create yk_dh_button
            app.yk_dh_button = uibutton(app.Tab_4, 'push');
            app.yk_dh_button.ButtonPushedFcn = createCallbackFcn(app, @yk_dh_buttonButtonPushed, true);
            app.yk_dh_button.Position = [108 20 100 26];
            app.yk_dh_button.Text = {'生成灰度图动画'; ''};

            % Create Label_4
            app.Label_4 = uilabel(app.Tab_4);
            app.Label_4.HorizontalAlignment = 'right';
            app.Label_4.Position = [48 60 77 22];
            app.Label_4.Text = '变化的物理量';

            % Create yk_DropDown
            app.yk_DropDown = uidropdown(app.Tab_4);
            app.yk_DropDown.Items = {'波长', '圆孔半径', ''};
            app.yk_DropDown.Position = [140 60 100 22];
            app.yk_DropDown.Value = '波长';

            % Create yk_1d_button
            app.yk_1d_button = uibutton(app.Tab_4, 'push');
            app.yk_1d_button.ButtonPushedFcn = createCallbackFcn(app, @yk_1d_buttonButtonPushed, true);
            app.yk_1d_button.Position = [20 94 75 26];
            app.yk_1d_button.Text = '强度曲线图';

            % Create yk_3d_button
            app.yk_3d_button = uibutton(app.Tab_4, 'push');
            app.yk_3d_button.ButtonPushedFcn = createCallbackFcn(app, @yk_3d_buttonButtonPushed, true);
            app.yk_3d_button.Position = [220 94 75 26];
            app.yk_3d_button.Text = '三维示意图';

            % Create Tab_5
            app.Tab_5 = uitab(app.TabGroup);
            app.Tab_5.Title = '设置';

            % Create Label_5
            app.Label_5 = uilabel(app.Tab_5);
            app.Label_5.Position = [42 452 65 22];
            app.Label_5.Text = '图样解析度';

            % Create Label_6
            app.Label_6 = uilabel(app.Tab_5);
            app.Label_6.Position = [142 452 53 22];
            app.Label_6.Text = '动画系数';

            % Create set_p_slider
            app.set_p_slider = uislider(app.Tab_5);
            app.set_p_slider.Limits = [100 2000];
            app.set_p_slider.Orientation = 'vertical';
            app.set_p_slider.ValueChangedFcn = createCallbackFcn(app, @set_p_sliderValueChanged, true);
            app.set_p_slider.Position = [42 136 3 277];
            app.set_p_slider.Value = 500;

            % Create set_a_slider
            app.set_a_slider = uislider(app.Tab_5);
            app.set_a_slider.Limits = [5 50];
            app.set_a_slider.Orientation = 'vertical';
            app.set_a_slider.ValueChangedFcn = createCallbackFcn(app, @set_a_sliderValueChanged, true);
            app.set_a_slider.Position = [142 136 3 277];
            app.set_a_slider.Value = 20;

            % Create set_p_spin
            app.set_p_spin = uispinner(app.Tab_5);
            app.set_p_spin.Limits = [100 2000];
            app.set_p_spin.RoundFractionalValues = 'on';
            app.set_p_spin.ValueDisplayFormat = '%d';
            app.set_p_spin.ValueChangedFcn = createCallbackFcn(app, @set_p_spinValueChanged, true);
            app.set_p_spin.HorizontalAlignment = 'center';
            app.set_p_spin.Position = [42 427 65 22];
            app.set_p_spin.Value = 500;

            % Create set_a_spin
            app.set_a_spin = uispinner(app.Tab_5);
            app.set_a_spin.Limits = [5 50];
            app.set_a_spin.RoundFractionalValues = 'on';
            app.set_a_spin.ValueDisplayFormat = '%d';
            app.set_a_spin.ValueChangedFcn = createCallbackFcn(app, @set_a_spinValueChanged, true);
            app.set_a_spin.HorizontalAlignment = 'center';
            app.set_a_spin.Position = [142 427 51 22];
            app.set_a_spin.Value = 20;

            % Create set_button
            app.set_button = uibutton(app.Tab_5, 'push');
            app.set_button.ButtonPushedFcn = createCallbackFcn(app, @set_buttonButtonPushed, true);
            app.set_button.Position = [108 11 100 26];
            app.set_button.Text = '保存';

            % Create set_Label
            app.set_Label = uilabel(app.Tab_5);
            app.set_Label.VerticalAlignment = 'top';
            app.set_Label.Position = [18 35 309 58];
            app.set_Label.Text = {'说明：'; '图样解析度越大，图样越清晰，过大可能导致运行缓慢；'; '动画系数越大，变量切换越慢，过大可能导致运行缓慢。'};

            % Create set_defult_Button
            app.set_defult_Button = uibutton(app.Tab_5, 'push');
            app.set_defult_Button.ButtonPushedFcn = createCallbackFcn(app, @set_defult_ButtonPushed, true);
            app.set_defult_Button.Position = [108 94 100 26];
            app.set_defult_Button.Text = '重置为默认值';

            % Create Label_7
            app.Label_7 = uilabel(app.Tab_5);
            app.Label_7.HorizontalAlignment = 'center';
            app.Label_7.Position = [233 258 65 22];
            app.Label_7.Text = '坐标轴显示';

            % Create axis_switch
            app.axis_switch = uiswitch(app.Tab_5, 'slider');
            app.axis_switch.Items = {'开', '关'};
            app.axis_switch.ValueChangedFcn = createCallbackFcn(app, @axis_switchValueChanged, true);
            app.axis_switch.Position = [242 283 45 20];
            app.axis_switch.Value = '开';

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;

            % Create UIAxes
            app.UIAxes = uiaxes(app.RightPanel);
            title(app.UIAxes, '夫琅禾费衍射仿真图样')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            app.UIAxes.TitleFontWeight = 'bold';
            app.UIAxes.Position = [1 6 542 549];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = myapp

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end