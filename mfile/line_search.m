function step = line_search(old,update,Data,para)

step_start = para.step_size(end)*1.3;%magic number
%step_start = 2;
%step_start = para.step_size(1);
tau = 0.8;
max_try = 15;
step = step_start;

cost_old = para.Cost.totalCost(end);

for i=1:max_try
    
    new = old + step*update;
    fidelity_new = compute_fidelity_for_line_search_yt(new,Data,para);
    cost_new = Cost_STCR(fidelity_new,new,para);

    if cost_new > cost_old
        step = step*tau;
    else
        %fprintf(['Step = ' num2str(step) '...\n'])
        %fprintf(['Cost = ' num2str(round(cost_new)) '...\n'])
        return
    end

end
%fprintf(['Step = ' num2str(step) '...\n'])
%fprintf(['Cost = ' num2str(round(cost_new)) '...\n'])