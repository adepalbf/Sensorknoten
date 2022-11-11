 classdef array
    
    properties
        dx
        dy
        pos
        
        phi_T
        Nx
        Ny
%         b_no_noise
%         SNR = 10
        ULA
    end
    
    properties(Dependent)
        A
        T
        principal_axis
    end
    methods
        function obj = array (varargin)
            % Definition of array microphone positions
            % center position and node orientation in degree off north
            %Parser Options
            p=inputParser;
            p.FunctionName='array';
            addOptional(p, 'pos',[0,0,0],@(x)validateattributes(x, {'numeric'},{'nonempty','real','vector'}));
            addOptional(p, 'phi_T',0,@(x)validateattributes(x, {'numeric'},{'nonempty','real','scalar'}));
            addOptional(p, 'dx',0.15,@(x)validateattributes(x, {'numeric'},{'nonempty','real','scalar','nonnegative'}));
            addOptional(p, 'dy',0.15,@(x)validateattributes(x, {'numeric'},{'nonempty','real','scalar','nonnegative'}));
            addOptional(p, 'Nx',2,@(x)validateattributes(x, {'numeric'},{'nonempty','real','scalar','integer','positive'}));
            addOptional(p, 'Ny',1,@(x)validateattributes(x, {'numeric'},{'nonempty','real','scalar','integer','positive'}));
            
            parse(p, varargin{:});
            
            obj.Nx = p.Results.Nx;
            obj.Ny = p.Results.Ny;
            obj.dx = p.Results.dx;
            obj.dy = p.Results.dy;
            obj.pos = p.Results.pos;
            obj.phi_T = p.Results.phi_T;
            
            N = obj.Nx * obj.Ny;
            obj.ULA = zeros(3,N);
            for ii = 1:obj.Nx
                for jj = 1:obj.Ny
                    obj.ULA(1,jj+obj.Ny*(ii-1)) = (ii-1)*obj.dx;
                    obj.ULA(2,jj+obj.Ny*(ii-1)) = (jj-1)*obj.dy;
                end
            end
            
        end
        
        function principal_axis = get.principal_axis(obj)
            principal_axis = obj.T*[1;0;0];
        end
        function T = get.T(obj)
            T = [cos(obj.phi_T), -sin(obj.phi_T), 0; sin(obj.phi_T), cos(obj.phi_T) 0; 0 0 1];
        end
        function A = get.A(obj)
            A = obj.T*obj.ULA;
            A = A+obj.pos';
        end
        
        function dn_max = dn_max(obj, fs)
            v = 343;
            
            difference = obj.A(:,2:end)-obj.A(:,1);
            distance = diag(sqrt(difference'*difference));
            
            % maximum TDOA in samples
            dn_max = round(fs*distance/v);
        end
        
        function phi_resolve = phi_resolve (obj, fs)
            phi_resolve = asin(1./obj.dn_max(fs))*180/pi;
        end
        
        function y = record(obj, fs,x,pos_source,amplitude_source, SNR, b_no_noise)
            y = nktp_sim(fs,x,pos_source,amplitude_source,SNR,b_no_noise, obj.A);
        end
        
        function [phi, P] = estimate_AOA(obj,y,fs, method)
            switch method
                case 'AOA_CC'
                    [beta, Detected_powers] = AOA_CC(y, obj.dn_max(fs));
                case 'AOA_GCC'
                    [beta, Detected_powers] = AOA_GCC(y, obj.dn_max(fs));
%                 case 'SAMV0'
%                     % steering vector matrix w.r.t all possible scanning DOA's
%                     DOAscan = -90: 0.5 :90; %0: 0.2: 180; % all possible DOA angles
%                     k = pi/obj.dx*[cos(DOAscan*pi/180); sin(DOAscan*pi/180)]; % sources in columns
%                     steering_array = exp(1i*mics*k);
%                     t_samples = 120;
%                     M = obj.Nx * obj.Ny;
%                     modulus_hat_das  = sum(abs(steering_array'*y/M), 2)/t_samples;
%                     [beta, Detected_powers] = SAMV0(y, steering_array, modulus_hat_das);
                case 'AOA_MUSIC'
                    [beta, Detected_powers] = AOA_MUSIC(y, obj.ULA);
                case 'AOA_ESPRITE'
                    [beta, Detected_powers] = AOA_ESPRITE(y, obj.ULA);
                case 'AOA_GCC_PHAT'
                    [beta, Detected_powers] = AOA_GCC_PHAT(y, obj.dn_max(fs));
            end
            
            phi = correct_angle(obj, beta);
            P = mean(Detected_powers);
        end
        
        function beta_av = correct_angle(obj, beta)
            
            
%             r_source = 5;
            
%             pos_source = r_source*[cos(phi_source), sin(phi_source),0];
%             plot(pos_source(1),pos_source(2), 'gx', 'LineWidth', 2)
            
            % raw estimate
%             beta_pos = obj.global_coords(beta, r_source);
%             plot([obj.pos(1)+zeros(size(beta_pos,1)), beta_pos(:,1)]',[obj.pos(2)+zeros(size(beta_pos,1)), beta_pos(:,2)]', ':')
                        
            difference = obj.A(:,2:end)-obj.A(:,1);
            
            angle_correction_factor = acos((difference'*obj.principal_axis)./sqrt(sum(difference.^2)'));

            % corrected estimate
            phi1 = mod(beta + angle_correction_factor, 2*pi);
            phi2 = mod(-beta + angle_correction_factor, 2*pi);
            
%             phi1_pos = obj.global_coords(phi1, r_source);
%             phi2_pos = obj.global_coords(phi2, r_source);
%             plot([obj.pos(1)+zeros(size(phi1_pos,1)), phi1_pos(:,1)]',[obj.pos(2)+zeros(size(phi1_pos,1)), phi1_pos(:,2)]', '--')
%             plot([obj.pos(1)+zeros(size(phi2_pos,1)), phi2_pos(:,1)]',[obj.pos(2)+zeros(size(phi2_pos,1)), phi2_pos(:,2)]', '--')
            
            C = idxcomb(phi1, phi2);
            
            D = abs(C(:,1)-C(:,2))+abs(C(:,1)-C(:,3));
            idx = find(D==min(D), 1);
            beta_av = mean(C(idx, :));
%             beta_av_pos = obj.global_coords(beta_av, r_source);
            
%             plot([obj.pos(1), beta_av_pos(1)],[obj.pos(2), beta_av_pos(2)], 'r', 'LineWidth', 2)
            
            
            %             xlim([0, r_source])
            %             ylim([-r_source, r_source])
%             beta_out = beta_av+obj.phi_T;
            
            %fprintf('beta_av = %.2f\n', beta_av*180/pi)
%             fprintf('beta_out = %.2f\n', beta_out*180/pi)
            %fprintf('beta_soll = %.2f\n', beta_soll*180/pi)
            
        end
        
        function error = error(obj, beta_est, phi_source, r_source)
            pos_source = r_source*[cos(phi_source), sin(phi_source),0];
            difference = pos_source-obj.pos;
            beta_soll = acos((difference*obj.principal_axis)./sqrt(sum(difference.^2)));
            error = abs(beta_soll-beta_est);
        end
        
        function coords = global_coords(obj, phi, r)
            coords = r*[cos(phi),sin(phi), zeros(size(phi))];
            for ii = 1:length(phi)
                coords(ii,:) = obj.T*coords(ii,:)';
                coords(ii,:) = coords(ii,:) + obj.pos;
            end
        end
        
        function plot_array(obj)            
%             C = hsv(size(obj.A, 2));
            
%             scatter(obj.A(1,:), obj.A(2,:), [], C)
            scatter(obj.A(1,:), obj.A(2,:))
            hold on
%             plot(obj.pos(1)+[0;obj.principal_axis(1)],obj.pos(2)+[0;obj.principal_axis(2)])
            axis equal
        end
    end
end