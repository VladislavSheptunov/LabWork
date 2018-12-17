package workshop.client.common;

import java.awt.Color;
import java.awt.Component;
import java.awt.Point;
import java.util.Arrays;
import javax.swing.JTable;
import javax.swing.table.DefaultTableCellRenderer;

public class SheduleTableRenderer extends DefaultTableCellRenderer {
    int size;
    Point[] filledCells;
    
    public SheduleTableRenderer(int size, Point[] filledCells) {
        this.size = size;
        this.filledCells = new Point[size];
        System.arraycopy(filledCells, 0, this.filledCells, 0, size);
    }

    @Override
    public Component getTableCellRendererComponent(JTable table, Object value,
            boolean isSelected, boolean hasFocus, int row, int column) {
        if (table == null) {
            return this;
        }
        if(Arrays.asList(filledCells).contains(new Point(row, column)))
            super.getTableCellRendererComponent(table, value, isSelected, 
                    hasFocus, row, column).setBackground(Color.red);
        else
            super.getTableCellRendererComponent(table, value, isSelected, 
                    hasFocus, row, column).setBackground(Color.white);
        return this;
    }
}
