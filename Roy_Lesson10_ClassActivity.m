% Define T_sd (Desired pose of the matrix provided)
T_sd = [-0.7071 -0.7071 0 0.7071;
        0.7071 -0.7071 0 2.1213;
        0 0 1 0;
        0 0 0 1];

% Define Lengths
L1 = 2;
L2 = 1;
t1 = 0;
t2 = 0;

%Definen Maximum iteration
max_iteration=1000;
tolerance=0.0001;

% Define q
q = [t1; t2];

for i=1:max_iteration

    % Update t1 and t2 from the latest q
    t1 = q(1);
    t2 = q(2);

    % Define T_sb
    T_sb = [cos(t1 + t2) -sin(t1 + t2) 0 L2*cos(t1+t2)+L1*cos(t1);
            sin(t1 + t2) cos(t1 + t2) 0 L2*sin(t1+t2)+L1*sin(t1);
            0 0 1 0;
            0 0 0 1];
    
    % Define the inverse of T_sb
    T_inv_sb = inv(T_sb);
    
    % Define R and its inverse
    R = T_sb(1:3, 1:3);
    R_inv = inv(R);
    
    % Extract the translation vector
    p = T_sb(1:3, 4);
    
    % Define T_bs as specified
    T_bs = [R_inv, -R_inv*p;
            0 0 0 1];
    
    % Define the skew-symmetric matrix [p]_x using the translation vector p
    p_x = [0, -p(3), p(2); 
           p(3), 0, -p(1); 
           -p(2), p(1), 0];
    
    % Construct the adjoint matrix Ad_T_bs
    Ad_T_bs = [R_inv, zeros(3,3);           
               p_x*R_inv, R_inv];
    
    %JACOBIAN CALCULATION
    % Obtain S1 and S2. Sw_i and a_i are obtained by visual inspection of the
    % robot picture in the lesson
    Sw1 = [0; 0; 1];
    Sw2 = [0; 0; 1];
    a1 = [0 0 0];
    a2 = [L1 0 0];
    Sv1 = transpose(cross(a1, Sw1));
    Sv2 = transpose(cross(a2, Sw2));
    S1 = [Sw1; Sv1];
    S2 = [Sw2; Sv2];
    
    % Obtain the T's for each joint. First obtain the Sw_brackets
    % (skew-symmetric representations of Sw_i)
    Sw1_bracket = [0 -Sw1(3) Sw1(2); 
                   Sw1(3) 0 -Sw1(1);
                   -Sw1(2) Sw1(1) 0];
    Sw2_bracket = [0 -Sw2(3) Sw2(2); 
                   Sw2(3) 0 -Sw2(1); 
                   -Sw2(2) Sw2(1) 0];
    
    % Obtain the R and p for each joint with Rodrigues' formula
    % JOINT 1
    R1 = eye(3) + sin(t1)*Sw1_bracket + (1-cos(t1))*(Sw1_bracket*Sw1_bracket); 
    p1 = (eye(3)*t1+(1-cos(t1))*Sw1_bracket+(t1-sin(t1))*(Sw1_bracket*Sw1_bracket))*Sv1;
    
    % JOINT 2
    R2 = eye(3) + sin(t2)*Sw2_bracket + (1-cos(t2))*(Sw2_bracket*Sw2_bracket);
    p2 = (eye(3)*t2+(1-cos(t2))*Sw2_bracket+(t2-sin(t2))*(Sw2_bracket*Sw2_bracket))*Sv2;
    
    
    % Obtain the Space Jacobian of each joint
    Js1 = S1; % Because it's the first joint
    p1_bracket = [0 -p1(3) p1(2); 
                  p1(3) 0 -p1(1); 
                  -p1(2) p1(1) 0]; % To obtain the Adjoint of T1
    
    zero_matrix = zeros(3); % To obtain the adjoint of T1
    
    AdT1 = [R1 zero_matrix; 
            p1_bracket*R1 R1];
    Js2 = AdT1*S2;
    
    % The final space Jacobian for this robot is
    Js = [Js1 Js2];
    
    %Define Body Jacobian from Adjoint of Tbs and Space Jacobian 
    J_b = Ad_T_bs * Js;
    
    %Define pseudo-inverse of Body Jacobian
    J_b_pseudo_inverse = pinv(J_b);
    
    % Define twist bracket, Nu bracket
    Vb_bracket = logm(T_inv_sb * T_sd); 
    
    %Extract Angular Velocity,from Nu bracket
    wb_bracket = Vb_bracket(1:3, 1:3); 
    wb_x = wb_bracket(3, 2);
    wb_y = wb_bracket(1, 3);
    wb_z = wb_bracket(2, 1);
    wb = [wb_x; wb_y; wb_z];
    
    %Extract Linear Velocity from Nu bracket
    vb = Vb_bracket(1:3, 4);
    Vb = [wb; vb];
    
    %Update the joint angles
    q_new = q+(J_b_pseudo_inverse*Vb);

    %Define error 
    error = norm(q_new - q);

    if error < tolerance || i == max_iteration
        break;
    end
    % Update q for the next iteration
    q = q_new;
end

% Final joint angles after convergence or reaching maximum iterations
disp('Final joint angles:');
disp(q);

