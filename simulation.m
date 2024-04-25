%Francis Mutetwa
%Student No: MTTFRA005
%Date: 21/03/2024
%Revised: 24/04/2024
%The following program simulates TOPSIS RAT selection. It does this by first 
%prompting the user for the number of users in a particular
%heterogenous network
%The heterogeneous network is already set to different RATs with certain
%criteria such as Available Bandwidth, Cost per Byte and Packet Delay.
% Prompt for number of calls in heterogeneous network
users = input('Enter the number of users in Network: ');
RAT1= 0;
RAT2= 0;
RAT3= 0;
for user = 1:users
    Selection_Criteria = {'Available Bandwidth', 'Cost(per Byte)', 'Delay(ms)'};

    % RATs data (normalized decision matrix)
    RAT_data = [
        0.969 0.628 0.510;  % 4G (RAT 1)
        0.242 0.189 0.845;  % Wi-Fi (RAT 2)
        0.032 0.754 0.169   % 3G (RAT 3)
    ];

    %User weights defined for each criterion
    weights = [randi([1,6]) randi([1,10]) randi([7,10])]; %random weights assigned
    %weights = [3,8,10]; % user weights more important for 3rd criterion

    %Normalising the weight matrix
    s = sum(weights);
    norm_weights = weights / s; %Nomalised weights

    %Calculating the weight normalised decision matrix
    w_n = RAT_data .*  norm_weights;

    %Determine the benefit and cost criteria
    %Determine the ideal solution (best), A* (best) and the 
    %negative ideal solution A- (worst) of H
    ideal_best = zeros(1, size(w_n, 2));
    ideal_worst = zeros(1, size(w_n, 2));

    for i = 1:3
        column = w_n(:,i);
        %cost criteria
        if strcmp(Selection_Criteria{i}, 'Cost(per Byte)') || strcmp(Selection_Criteria{i}, 'Delay(ms)')
            ideal_best(i)=min(column);
            ideal_worst(i)=max(column);

        %benefit criteria
        elseif strcmp(Selection_Criteria{i}, 'Available Bandwidth') 
            ideal_best(i)=max(column);
            ideal_worst(i)=min(column);
        end
    end

    % Initialize an array to store the positive ideal separation(d*)
    positive_ideal_separation  = zeros(3,1);
    negative_ideal_separation  = zeros(3,1);
    % Compute the squared differences for each row of matrix x
    for i = 1:3
        positive_ideal_separation(i) = sqrt(sum((w_n(i,:) - ideal_best).^2));
        negative_ideal_separation(i) = sqrt(sum((w_n(i,:) - ideal_worst).^2));
    end

    %Calculating the closeness coefficient of each of the available 
    %RATs to the ideal solution, and negative ideal solution.

    % Initialize an array to store the closeness coefficients
    closeness_coefficients = zeros(3,1);
    for n = 1:3
        closeness_coefficients(n) =  negative_ideal_separation(n) / ( positive_ideal_separation(n) + negative_ideal_separation(n));
    end
    %Finding the RAT with the highest value
    rat = find(closeness_coefficients == max(closeness_coefficients));

    %Calculating the distribution of calls
    if rat==1
        RAT1 = RAT1+1;
    elseif rat ==2
        RAT2 = RAT2+1;
    else
        RAT3 = RAT3+1;
    end
end


%Show results on bar graph.
call_distro=[RAT1, RAT2, RAT3];
bar(call_distro);
ylabel("Number of Calls");
RAT_labels = {'4G', 'Wi-Fi', '3G'};
xticklabels(RAT_labels);
title("Distribution of Calls");





