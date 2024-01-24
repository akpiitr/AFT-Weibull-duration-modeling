function [prob_vals] = deterioration_smp(CR_i, planning_horizon, AADTs_i)

CR=CR_i;
T = planning_horizon;
dt = 1;
time_steps = 1:dt:T;

%% calculating prob CR->CR
if CR>4
    p_9to9 = [];
    for t=1:length(time_steps)
        [pdf_9, cdf_9, haz_9] = pdf_cdf_haz(time_steps(t), CR, AADTs_i(time_steps(t)));
        p_9to9(t) = 1-cdf_9;
    end
end
%% calculating prob CR->CR-1   
if CR>4
    p_9to7a =[];
    ex = [];
    for t=1:length(time_steps)
        dls=[];
        for t2=1:dt:time_steps(t)
            [pdf_9, cdf_9, haz_9] = pdf_cdf_haz(t2, CR, AADTs_i(t2));
            [pdf_8, cdf_8, haz_8] = pdf_cdf_haz(time_steps(t)-t2, CR-1, AADTs_i(t));
            val = pdf_9*cdf_8;
            dls = [dls; val];
            ex = [ex;val];
        end
        p_9to7a = [p_9to7a; sum(dls,'all')];
    end
    p_9to7a = p_9to7a';
    p_9to8 = 1-p_9to9-p_9to7a;
end
%% calculating CR->CR-2
if CR>5
    p_9to6a =[];
    for t=1:length(time_steps)
        val = 0;
        for a = 1:t
            for b = a+1:t
                if  a~=b && b~=t
                    [pdf_9, cdf_9, haz_9] = pdf_cdf_haz(a, CR, AADTs_i(a));
                    [pdf_8, cdf_8, haz_8] = pdf_cdf_haz(b-a, CR-1, AADTs_i(b-a));
                    [pdf_7, cdf_7, haz_7] = pdf_cdf_haz(t-b, CR-2, AADTs_i(t-b));
                    val = val + pdf_9 * pdf_8 * cdf_7;
                end
            end
        end
        p_9to6a = [p_9to6a; val];
    end
    p_9to6a = p_9to6a';
    p_9to7 = 1-p_9to9-p_9to8-p_9to6a;
end   
%% calculating CR->CR-3
if CR>6
    p_9to5a =[];
    for t=1:length(time_steps)
        val = 0;
        for a = 1:t
            for b = a+1:t
                for c= b:t
                    if  a~=b && b~=c && t~=c
                        [pdf_9, cdf_9, haz_9] = pdf_cdf_haz(a, CR, AADTs_i(a));
                        [pdf_8, cdf_8, haz_8] = pdf_cdf_haz(b-a, CR-1, AADTs_i(b-a));
                        [pdf_7, cdf_7, haz_7] = pdf_cdf_haz(c-b, CR-2, AADTs_i(c-b));
                        [pdf_6, cdf_6, haz_6] = pdf_cdf_haz(t-c, CR-3, AADTs_i(t-c));
                        val = val + pdf_9 * pdf_8* pdf_7 * cdf_6;
                    end
                end
            end
        end
        p_9to5a = [p_9to5a; val];
    end
    p_9to5a = p_9to5a';
    p_9to6 = 1- p_9to9 - p_9to8 - p_9to7 - p_9to5a;
end    
%% calculating CR->CR-4
if CR>7
    p_9to4a =[];
    for t=1:length(time_steps)
        val = 0;
        for a = 1:t
            for b = a:t
                for c= b:t
                    for d= c:t
                        if  a~=b && b~=c && d~=c && t~=d
                            [pdf_9, cdf_9, haz_9] = pdf_cdf_haz(a, CR, AADTs_i(a));
                            [pdf_8, cdf_8, haz_8] = pdf_cdf_haz(b-a, CR-1, AADTs_i(b-a));
                            [pdf_7, cdf_7, haz_7] = pdf_cdf_haz(c-b, CR-2, AADTs_i(c-b));
                            [pdf_6, cdf_6, haz_6] = pdf_cdf_haz(d-c, CR-3, AADTs_i(d-c));
                            [pdf_5, cdf_5, haz_5] = pdf_cdf_haz(t-d, CR-4, AADTs_i(t-d));
                            val = val + pdf_9 * pdf_8* pdf_7 * pdf_6 * cdf_5;
                        end
                    end
                end
            end
        end
        p_9to4a = [p_9to4a; val];
    end
    p_9to4a = p_9to4a';
    p_9to5 = 1- p_9to9 - p_9to8 - p_9to7 - p_9to6 - p_9to4a;
end    
%%
prob_vals=[];  %(transition_prob, time)
if CR==9
    prob_vals = [p_9to9; p_9to8; p_9to7; p_9to6; p_9to5; p_9to4a];
elseif CR==8
    prob_vals = [p_9to9; p_9to8; p_9to7; p_9to6; p_9to5a];
elseif CR==7
    prob_vals = [p_9to9; p_9to8; p_9to7; p_9to6a];
elseif CR==6
    prob_vals = [p_9to9; p_9to8; p_9to7a];
elseif CR==5
    prob_vals = [p_9to9; 1-p_9to9];
elseif CR<=4
    prob_vals = [zeros(1,T); ones(1,T)];
end

end
