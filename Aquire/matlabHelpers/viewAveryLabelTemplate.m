function [] = viewAveryLabelTemplate(I,templateData)

    imshow(I,[]);
    hold on
    for e = 1:size(templateData.labelBox,1)
        rectangle('Position',templateData.labelBox(e,:),'EdgeColor','r');
        rectangle('Position',templateData.qrBox(e,:),'EdgeColor','g');
        rectangle('Position',templateData.textBox(e,:),'EdgeColor','b');
        text(templateData.centroid(e,1),templateData.centroid(e,2),num2str(e));
    end

end