% test_P2p1_Granberry
format compact;

x_list_from_GE = [];
x_list_from_rref = [];
sys_list = [];
for verification_matries = 1:5 % Change verification_matries = 1:6 to see non_square matrix error
    switch verification_matries
        case 1
            A = [3 4 1;
                2 -2 1;
                5 4 0.4];
            
            b = [7; -1; 9];
            
        case 2
            A = [0 1 1;
                1 -2 1;
                7 2 2];
            
            b = [12; 1; 8];
            
        case 3
            A = [1 -1 2;
                0 4 0;
                0 2 1];
            
            b = [22; 44; 9];
            
        case 4
            A = [0.0005 1;
                1 1];
            
            b = [1; 2];
            
        case 5
            A = [0 -1 2 1 1;
                1 -1 1 1 1;
                2 1 3 2 2;
                2 -3 -4 -3 0;
                -1 5 0 0 -1];
            
            b = [4; 4; 12; 7; -1];
            
        case 6 % A is not a square.
            A = [1 2 3;
                5 6 9];
            
            b = [4 1]';
    end
    
    % gauss_elim_P2p1
    x_from_GE = gauss_elim_P2p1(A,b);
    x_from_GE = x_from_GE';
    
    % rref
    aug_from_rref = rref([A b]);
    x_from_rref = aug_from_rref(:,length(aug_from_rref));
    x_from_rref = x_from_rref';
    
    % Concatenating x lists
    x_list_from_GE = [x_list_from_GE  x_from_GE];
    x_list_from_rref = [x_list_from_rref x_from_rref];
    
    % Creating a list of Sys#_x#
    for n = 1:length(x_from_rref)
        x = num2str(double(n));
        s = ['Sys' num2str(double(verification_matries)) '_x' x];
        sys_list = [sys_list; s];
    end
end
Vector_x_talbe = table(sys_list, x_list_from_GE', x_list_from_rref');
Vector_x_talbe.Properties.VariableNames = {' ', 'x_from_GE', 'x_from_rref(aug)'};
Vector_x_talbe