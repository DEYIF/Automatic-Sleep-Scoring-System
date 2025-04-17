function [precision, recall, f1, specificity] = calcMetrics(true_labels, pred_labels)
    % positive 1，negative 0
    TP = sum((pred_labels == 1) & (true_labels == 1));
    TN = sum((pred_labels == 0) & (true_labels == 0));
    FP = sum((pred_labels == 1) & (true_labels == 0));
    FN = sum((pred_labels == 0) & (true_labels == 1));
    
    precision = TP / (TP + FP);
    recall = TP / (TP + FN);
    f1 = 2 * precision * recall / (precision + recall);
    specificity = TN / (TN + FP);
end