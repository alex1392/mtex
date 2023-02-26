function plotPoleFigures(ori, CS, xlabel, ylabel)
figure
if strcmp(CS.pointGroup,'6/mmm')
	plotPDF(ori,...
		[Miller({0 0 0 1}, CS) Miller({1 1 -2 0},CS)],...
		'contour',100,...
		'contourf', 'linecolor','none','minmax')
elseif strcmp(CS.pointGroup, 'm-3m')
	plotPDF(ori,...
		[Miller({1 1 0},CS) Miller({1 1 1},CS)],...
		'contour',100,...
		'contourf', 'linecolor','none','minmax')
end
setAxis();

function setAxis()
fig = gcf;
for i = 1:length(fig.Children)
	if isa(fig.Children(i),'matlab.graphics.axis.Axes')
		ax = fig.Children(i);
		text(vector3d.X, xlabel,'horizontalAlignment','left', 'parent', ax)
		text(vector3d.Y, ylabel,'VerticalAlignment','bottom', 'parent', ax)
		ax.Title.Position = [0 -1.7 0];
	end
end
setColorRange('equal')
mtexColorbar
mtexColorMap BlueJet
end
end