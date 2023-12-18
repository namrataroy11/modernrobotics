% Clean the console
clear; clc

% Define our symbolic variables
syms t1(t) t2(t) L1 L2

% Define q and q_dot
q = [t1; t2];
q_dot = simplify(diff(q, t));

% Define T (the T for this robot is in the wiki for Lesson 9)
T = [cos(t1 + t2) -sin(t1 + t2) 0 L2*cos(t1+t2)+L1*cos(t1);
    sin(t1 + t2) cos(t1 + t2) 0 L2*sin(t1+t2)+L1*sin(t1);
    0 0 1 0;
    0 0 0 1];

% Obtain the first derivative of T and the inverse of T (for use later)
T_dot = simplify(diff(T, t));
T_inv = simplify(inv(T));

% FIRST METHOD OF OBTAINING OF THE SPACE TWIST Vs
Vs_bracket = simplify(T_dot * T_inv); % This returns a symfun object, and it is not indexable
Vs_bracket = Vs_bracket(t); % This evaluates the symfun object, return a sym object, and that IS indexable!
ws_bracket = Vs_bracket(1:3, 1:3); % This indexing would not work without evaluating the symfun object first!
vs = Vs_bracket(1:3, 4);
ws_x = ws_bracket(3, 2);
ws_y = ws_bracket(1, 3);
ws_z = ws_bracket(2, 1);
ws = [ws_x; ws_y; ws_z];
Vs = [ws; vs]; % FINAL ANSWER for this method

% SECOND METHOD OF OBTAINING THE SPACE TWIST Vs

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
R1 = eye(3) + sin(t1)*Sw1_bracket + (1-cos(t1))*(Sw1_bracket*Sw1_bracket); %symfun object, not indexable
R1 = R1(t); % sym object, now indexable!
p1 = (eye(3)*t1+(1-cos(t1))*Sw1_bracket+(t1-sin(t1))*(Sw1_bracket*Sw1_bracket))*Sv1;
p1 = p1(t);

% JOINT 2
R2 = eye(3) + sin(t2)*Sw2_bracket + (1-cos(t2))*(Sw2_bracket*Sw2_bracket);
R2 = R2(t);
p2 = (eye(3)*t2+(1-cos(t2))*Sw2_bracket+(t2-sin(t2))*(Sw2_bracket*Sw2_bracket))*Sv2;
p2 = p2(t);

% T's for each joint
row = [0 0 0 1];
T1_no_row = [R1 p1];
T2_no_row = [R2 p2];
T1 = [T1_no_row; row];
T2 = [T2_no_row; row];

% Obtain the Space Jacobian of each joint
Js1 = S1; % Because it's the first joint
p1_bracket = [0 -p1(3) p1(2); p1(3) 0 -p1(1); -p1(2) p1(1) 0]; % To obtain the Adjoint of T1
zero_matrix = zeros(3); % To obtain the adjoint of T1
AdT1 = [R1 zero_matrix; p1_bracket*R1 R1];
Js2 = AdT1*S2;

% The final space Jacobian for this robot is
Js = [Js1 Js2];

% The final Vs, with the second method, is:
Vs2 = Js*q_dot;

% Now to print all the answers
Vs;
Vs2;
Js;

