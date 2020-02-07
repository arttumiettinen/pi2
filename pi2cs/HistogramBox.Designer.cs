namespace pi2cs
{
    partial class HistogramBox
    {
        /// <summary> 
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary> 
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Component Designer generated code

        /// <summary> 
        /// Required method for Designer support - do not modify 
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            System.Windows.Forms.DataVisualization.Charting.ChartArea chartArea1 = new System.Windows.Forms.DataVisualization.Charting.ChartArea();
            System.Windows.Forms.DataVisualization.Charting.Series series1 = new System.Windows.Forms.DataVisualization.Charting.Series();
            this.chartHist = new System.Windows.Forms.DataVisualization.Charting.Chart();
            this.panel1 = new System.Windows.Forms.Panel();
            this.tableLayoutPanel1 = new System.Windows.Forms.TableLayoutPanel();
            this.label2 = new System.Windows.Forms.Label();
            this.label1 = new System.Windows.Forms.Label();
            this.scrollMaximum = new System.Windows.Forms.HScrollBar();
            this.scrollMinimum = new System.Windows.Forms.HScrollBar();
            this.flowLayoutPanel1 = new System.Windows.Forms.FlowLayoutPanel();
            this.checkFullRange = new System.Windows.Forms.CheckBox();
            this.labelMinimum = new System.Windows.Forms.Label();
            this.labelMaximum = new System.Windows.Forms.Label();
            ((System.ComponentModel.ISupportInitialize)(this.chartHist)).BeginInit();
            this.panel1.SuspendLayout();
            this.tableLayoutPanel1.SuspendLayout();
            this.flowLayoutPanel1.SuspendLayout();
            this.SuspendLayout();
            // 
            // chartHist
            // 
            chartArea1.AxisX.LabelStyle.Enabled = false;
            chartArea1.AxisY.LabelStyle.Enabled = false;
            chartArea1.Name = "ChartArea1";
            this.chartHist.ChartAreas.Add(chartArea1);
            this.chartHist.Dock = System.Windows.Forms.DockStyle.Fill;
            this.chartHist.Location = new System.Drawing.Point(0, 0);
            this.chartHist.Margin = new System.Windows.Forms.Padding(2);
            this.chartHist.Name = "chartHist";
            series1.ChartArea = "ChartArea1";
            series1.Name = "Series1";
            this.chartHist.Series.Add(series1);
            this.chartHist.Size = new System.Drawing.Size(377, 151);
            this.chartHist.TabIndex = 2;
            // 
            // panel1
            // 
            this.panel1.Controls.Add(this.tableLayoutPanel1);
            this.panel1.Dock = System.Windows.Forms.DockStyle.Bottom;
            this.panel1.Location = new System.Drawing.Point(0, 151);
            this.panel1.Margin = new System.Windows.Forms.Padding(2);
            this.panel1.Name = "panel1";
            this.panel1.Size = new System.Drawing.Size(377, 81);
            this.panel1.TabIndex = 3;
            // 
            // tableLayoutPanel1
            // 
            this.tableLayoutPanel1.ColumnCount = 3;
            this.tableLayoutPanel1.ColumnStyles.Add(new System.Windows.Forms.ColumnStyle());
            this.tableLayoutPanel1.ColumnStyles.Add(new System.Windows.Forms.ColumnStyle(System.Windows.Forms.SizeType.Percent, 100F));
            this.tableLayoutPanel1.ColumnStyles.Add(new System.Windows.Forms.ColumnStyle());
            this.tableLayoutPanel1.Controls.Add(this.label2, 0, 1);
            this.tableLayoutPanel1.Controls.Add(this.label1, 0, 0);
            this.tableLayoutPanel1.Controls.Add(this.scrollMaximum, 1, 1);
            this.tableLayoutPanel1.Controls.Add(this.scrollMinimum, 1, 0);
            this.tableLayoutPanel1.Controls.Add(this.flowLayoutPanel1, 1, 2);
            this.tableLayoutPanel1.Controls.Add(this.labelMinimum, 2, 0);
            this.tableLayoutPanel1.Controls.Add(this.labelMaximum, 2, 1);
            this.tableLayoutPanel1.Dock = System.Windows.Forms.DockStyle.Fill;
            this.tableLayoutPanel1.Location = new System.Drawing.Point(0, 0);
            this.tableLayoutPanel1.Margin = new System.Windows.Forms.Padding(2);
            this.tableLayoutPanel1.Name = "tableLayoutPanel1";
            this.tableLayoutPanel1.Padding = new System.Windows.Forms.Padding(4, 8, 8, 4);
            this.tableLayoutPanel1.RowCount = 3;
            this.tableLayoutPanel1.RowStyles.Add(new System.Windows.Forms.RowStyle());
            this.tableLayoutPanel1.RowStyles.Add(new System.Windows.Forms.RowStyle());
            this.tableLayoutPanel1.RowStyles.Add(new System.Windows.Forms.RowStyle(System.Windows.Forms.SizeType.Absolute, 16F));
            this.tableLayoutPanel1.Size = new System.Drawing.Size(377, 81);
            this.tableLayoutPanel1.TabIndex = 1;
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(6, 29);
            this.label2.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(27, 13);
            this.label2.TabIndex = 4;
            this.label2.Text = "Max";
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(6, 8);
            this.label1.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(24, 13);
            this.label1.TabIndex = 3;
            this.label1.Text = "Min";
            // 
            // scrollMaximum
            // 
            this.scrollMaximum.Dock = System.Windows.Forms.DockStyle.Top;
            this.scrollMaximum.LargeChange = 20;
            this.scrollMaximum.Location = new System.Drawing.Point(35, 29);
            this.scrollMaximum.Maximum = 400;
            this.scrollMaximum.Name = "scrollMaximum";
            this.scrollMaximum.Size = new System.Drawing.Size(320, 21);
            this.scrollMaximum.TabIndex = 2;
            this.scrollMaximum.ValueChanged += new System.EventHandler(this.scrollMaximum_ValueChanged);
            // 
            // scrollMinimum
            // 
            this.scrollMinimum.Dock = System.Windows.Forms.DockStyle.Top;
            this.scrollMinimum.LargeChange = 20;
            this.scrollMinimum.Location = new System.Drawing.Point(35, 8);
            this.scrollMinimum.Maximum = 400;
            this.scrollMinimum.Name = "scrollMinimum";
            this.scrollMinimum.Size = new System.Drawing.Size(320, 21);
            this.scrollMinimum.TabIndex = 1;
            this.scrollMinimum.ValueChanged += new System.EventHandler(this.scrollMinimum_ValueChanged);
            // 
            // flowLayoutPanel1
            // 
            this.flowLayoutPanel1.Controls.Add(this.checkFullRange);
            this.flowLayoutPanel1.Dock = System.Windows.Forms.DockStyle.Fill;
            this.flowLayoutPanel1.Location = new System.Drawing.Point(37, 52);
            this.flowLayoutPanel1.Margin = new System.Windows.Forms.Padding(2);
            this.flowLayoutPanel1.Name = "flowLayoutPanel1";
            this.flowLayoutPanel1.Size = new System.Drawing.Size(316, 23);
            this.flowLayoutPanel1.TabIndex = 5;
            // 
            // checkFullRange
            // 
            this.checkFullRange.AutoSize = true;
            this.checkFullRange.Location = new System.Drawing.Point(2, 2);
            this.checkFullRange.Margin = new System.Windows.Forms.Padding(2);
            this.checkFullRange.Name = "checkFullRange";
            this.checkFullRange.Size = new System.Drawing.Size(72, 17);
            this.checkFullRange.TabIndex = 6;
            this.checkFullRange.Text = "Full range";
            this.checkFullRange.UseVisualStyleBackColor = true;
            this.checkFullRange.CheckedChanged += new System.EventHandler(this.checkFullRange_CheckedChanged);
            // 
            // labelMinimum
            // 
            this.labelMinimum.AutoSize = true;
            this.labelMinimum.Location = new System.Drawing.Point(357, 8);
            this.labelMinimum.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.labelMinimum.Name = "labelMinimum";
            this.labelMinimum.Size = new System.Drawing.Size(10, 13);
            this.labelMinimum.TabIndex = 6;
            this.labelMinimum.Text = "-";
            // 
            // labelMaximum
            // 
            this.labelMaximum.AutoSize = true;
            this.labelMaximum.Location = new System.Drawing.Point(357, 29);
            this.labelMaximum.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.labelMaximum.Name = "labelMaximum";
            this.labelMaximum.Size = new System.Drawing.Size(10, 13);
            this.labelMaximum.TabIndex = 7;
            this.labelMaximum.Text = "-";
            // 
            // HistogramBox
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.Controls.Add(this.chartHist);
            this.Controls.Add(this.panel1);
            this.Name = "HistogramBox";
            this.Size = new System.Drawing.Size(377, 232);
            ((System.ComponentModel.ISupportInitialize)(this.chartHist)).EndInit();
            this.panel1.ResumeLayout(false);
            this.tableLayoutPanel1.ResumeLayout(false);
            this.tableLayoutPanel1.PerformLayout();
            this.flowLayoutPanel1.ResumeLayout(false);
            this.flowLayoutPanel1.PerformLayout();
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.DataVisualization.Charting.Chart chartHist;
        private System.Windows.Forms.Panel panel1;
        private System.Windows.Forms.TableLayoutPanel tableLayoutPanel1;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.HScrollBar scrollMaximum;
        private System.Windows.Forms.HScrollBar scrollMinimum;
        private System.Windows.Forms.FlowLayoutPanel flowLayoutPanel1;
        private System.Windows.Forms.CheckBox checkFullRange;
        private System.Windows.Forms.Label labelMinimum;
        private System.Windows.Forms.Label labelMaximum;
    }
}
