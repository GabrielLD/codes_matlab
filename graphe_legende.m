function [hx,hy,htitre]=graphe_legende(xtitre,ytitre,titre,graphe_multiple)

size=36;
min=0.25;
size_titre=36;

hx=xlabel(xtitre,'fontsize',size,'fontname','times','interpreter','latex');
hy=ylabel(ytitre,'fontsize',size,'fontname','times','interpreter','latex','rotation',90,'units','normalized','position',[-0.08 0.5]);
htitre=title(titre,'fontsize',size_titre,'fontname','times','interpreter','latex');

set(gca,'Box','on');


if (graphe_multiple)
 %   set(gcf,'OuterPosition',[200 500 1000 400])
%    set(gcf,'ActivePositionProperty','OuterPosition')
else
%set(gca,'fontname','Times','fontsize',size-2)
set(gca,'fontname','Times','fontsize',14,'units','normalized','position',[min-0.02 min-0.02 0.95-min 0.85-min]);
%set(gcf,'Position',[200 800 400 400]);
end

%%%%%%%%%%%% Pour changer les labels sur chaque axe
set(gca,'xtick',[0 100 200 300 400 500],'fontname','arial','fontsize',20)
set(gca,'ytick',[0 5 10 15 20 25],'fontname','arial','fontsize',20)

%h = axes('Position',[0.6 0.6 0.9 0.9],'Visible','off');
%axes('Position',[.25 .1 .7 .8])

