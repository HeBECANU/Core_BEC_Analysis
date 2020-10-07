function [plt,f] = ci_plot(X,Y,I,varargin)
% a widget to plot lines with some surrounding interval such as standard
% error or confidence interval.
    p = inputParser;
%     defaultFontSize = 12; 
%     defaultFont = 'times'; 
    addParameter(p,'mask',ones(size(X)));
    addParameter(p,'LogScale',0);
    addParameter(p,'FaceAlpha',0.2);
    addParameter(p,'Ymin',nan);
    addParameter(p,'LineCol',[0,0,0]);
    addParameter(p,'AreaCol',0.5*[1,1,1]);
    addParameter(p,'LineWidth',2);    
    addParameter(p,'mode','Interval');
    parse(p,varargin{:});
    areacol = p.Results.AreaCol;
    fill_alpha = p.Results.FaceAlpha;
    linecol = p.Results.LineCol;
    plotmode = p.Results.mode;
    ymin = p.Results.Ymin;
    lw = p.Results.LineWidth;
    logscale = p.Results.LogScale;
    
    X = col_vec(X);

    
    if any(size(I)==1)
        Y_upper = col_vec(Y)+col_vec(I);
        Y_lower = col_vec(Y)-col_vec(I);
%     end
    elseif any(size(I) == 2) %upper and lower CI specified
        if size(I,1) == 2
            I = I';
        end
        if strcmp(plotmode,'Interval') % supplied vals are the CI widths
            Y_upper = col_vec(Y)+col_vec(I(:,2));
            Y_lower = col_vec(Y)-col_vec(I(:,1));
        elseif strcmp(plotmode,'Absolute')
            Y_upper = I(:,1);
            Y_lower = I(:,2);
        else
            error('Plot mode specification error')
        end
    else
        error('CI size specification error')
    end
    
    if logscale
        if isnan(ymin)
            ymin = 0.1*min(Y(Y>0));
        end
        Y_lower(Y_lower<ymin) = ymin;
    end
    
    
    plt=plot(X,Y,'Color',linecol,'LineWidth',lw);
    hold on
    f=fill([X;flipud(X)],[Y_upper;flipud(Y_lower)],areacol,'LineStyle','none','FaceAlpha',fill_alpha);
    uistack(f,'bottom')
%     plot(X,Y_upper,':','Color',areacol)
%     plot(X,Y_lower,':','Color',areacol)
    
    if logscale
        ylim([ymin,2*max(Y_upper)])
        set(gca,'Yscale','log')
    end
%     fill([col_vec(X);flipud(col_vec(X))],[(col_vec(Y)+col_vec(I));flipud((col_vec(Y)-col_vec(I)))],...
%                 [0.5,0.5,0.5],'FaceAlpha',0.5,'LineStyle','none')

end