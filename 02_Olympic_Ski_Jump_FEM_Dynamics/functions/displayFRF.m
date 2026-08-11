function displayFRF(a,om, X, t, newFigure)
    if newFigure    
        figure();
    end 
    sgtitle(t);
    subplot(2,1,1);
    semilogy(om/2/pi,abs(squeeze(X(a,:))));
    hold on;
    grid on;
    ylabel('Amplitude [m]');
    subplot(2,1,2);
    plot(om/2/pi,unwrap(angle(X(a,:))));
    hold on;
    grid on;
    ylabel('phase [rad]');
    xlabel('frequency [Hz]');
end