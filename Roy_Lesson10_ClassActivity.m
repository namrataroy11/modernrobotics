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
    
    % Define the inverse of T_bs
    T_bs= inv(T_sb);
    
    % Define R and its inverse
    R = T_bs(1:3, 1:3);
    % Extract the translation vector
    p = T_bs(1:3, 4);
    p_bracket = [0, -p(3), p(2); 
           p(3), 0, -p(1); 
           -p(2), p(1), 0];

    zero_matrix = zeros(3);
    
    % Define T_bs as specified
    T_bs = [R, -R_inv*p;
            0 0 0 1];
    
    % Construct the adjoint matrix Ad_T_bs
    Ad_T_bs = [R, zero_matrix;           
               p_bracket*R, R];
    
    %JACOBIAN CALCULATION
    %Define Space Jacobian 
    Js = [0, 0; 
          0, 0;
          1, 1;
          0,  L1*sin(t1);
          0, -L1*cos(t1);
          0, 0];
    
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

