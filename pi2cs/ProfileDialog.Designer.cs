namespace pi2cs
{
    partial class ProfileDialog
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

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            System.Windows.Forms.DataVisualization.Charting.ChartArea chartArea1 = new System.Windows.Forms.DataVisualization.Charting.ChartArea();
            System.Windows.Forms.DataVisualization.Charting.Series series1 = new System.Windows.Forms.DataVisualization.Charting.Series();
            this.chart = new System.Windows.Forms.DataVisualization.Charting.Chart();
            this.tableLayoutPanel1 = new System.Windows.Forms.TableLayoutPanel();
            this.labelProfileWidth = new System.Windows.Forms.Label();
            this.numericUpDownProfileWidth = new System.Windows.Forms.NumericUpDown();
            ((System.ComponentModel.ISupportInitialize)(this.chart)).BeginInit();
            this.tableLayoutPanel1.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.numericUpDownProfileWidth)).BeginInit();
            this.SuspendLayout();
            // 
            // chart
            // 
            chartArea1.AxisX.LabelStyle.Format = "{0:0}";
            chartArea1.AxisX.Title = "Position [pix]";
            chartArea1.AxisY.Title = "Gray value";
            chartArea1.Name = "ChartArea";
            this.chart.ChartAreas.Add(chartArea1);
            this.chart.Dock = System.Windows.Forms.DockStyle.Fill;
            this.chart.Location = new System.Drawing.Point(0, 0);
            this.chart.Name = "chart";
            series1.ChartArea = "ChartArea";
            series1.ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Line;
            series1.Color = System.Drawing.Color.Red;
            series1.Name = "Profile";
            this.chart.Series.Add(series1);
            this.chart.Size = new System.Drawing.Size(501, 245);
            this.chart.TabIndex = 0;
            this.chart.Text = "chart1";
            // 
            // tableLayoutPanel1
            // 
            this.tableLayoutPanel1.ColumnCount = 2;
            this.tableLayoutPanel1.ColumnStyles.Add(new System.Windows.Forms.ColumnStyle());
            this.tableLayoutPanel1.ColumnStyles.Add(new System.Windows.Forms.ColumnStyle(System.Windows.Forms.SizeType.Percent, 100F));
            this.tableLayoutPanel1.Controls.Add(this.labelProfileWidth, 0, 0);
            this.tableLayoutPanel1.Controls.Add(this.numericUpDownProfileWidth, 1, 0);
            this.tableLayoutPanel1.Dock = System.Windows.Forms.DockStyle.Bottom;
            this.tableLayoutPanel1.Location = new System.Drawing.Point(0, 245);
            this.tableLayoutPanel1.Name = "tableLayoutPanel1";
            this.tableLayoutPanel1.RowCount = 1;
            this.tableLayoutPanel1.RowStyles.Add(new System.Windows.Forms.RowStyle(System.Windows.Forms.SizeType.Percent, 100F));
            this.tableLayoutPanel1.RowStyles.Add(new System.Windows.Forms.RowStyle(System.Windows.Forms.SizeType.Absolute, 30F));
            this.tableLayoutPanel1.Size = new System.Drawing.Size(501, 30);
            this.tableLayoutPanel1.TabIndex = 1;
            // 
            // labelProfileWidth
            // 
            this.labelProfileWidth.AutoSize = true;
            this.labelProfileWidth.Dock = System.Windows.Forms.DockStyle.Fill;
            this.labelProfileWidth.Location = new System.Drawing.Point(3, 0);
            this.labelProfileWidth.Name = "labelProfileWidth";
            this.labelProfileWidth.Size = new System.Drawing.Size(67, 30);
            this.labelProfileWidth.TabIndex = 0;
            this.labelProfileWidth.Text = "Profile width:";
            this.labelProfileWidth.TextAlign = System.Drawing.ContentAlignment.MiddleCenter;
            // 
            // numericUpDownProfileWidth
            // 
            this.numericUpDownProfileWidth.Anchor = System.Windows.Forms.AnchorStyles.Left;
            this.numericUpDownProfileWidth.Increment = new decimal(new int[] {
            2,
            0,
            0,
            0});
            this.numericUpDownProfileWidth.Location = new System.Drawing.Point(76, 5);
            this.numericUpDownProfileWidth.Minimum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.numericUpDownProfileWidth.Name = "numericUpDownProfileWidth";
            this.numericUpDownProfileWidth.Size = new System.Drawing.Size(55, 20);
            this.numericUpDownProfileWidth.TabIndex = 1;
            this.numericUpDownProfileWidth.Value = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.numericUpDownProfileWidth.ValueChanged += new System.EventHandler(this.numericUpDownProfileWidth_ValueChanged);
            // 
            // ProfileDialog
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(501, 275);
            this.Controls.Add(this.chart);
            this.Controls.Add(this.tableLayoutPanel1);
            this.FormBorderStyle = System.Windows.Forms.FormBorderStyle.SizableToolWindow;
            this.MaximizeBox = false;
            this.MinimizeBox = false;
            this.Name = "ProfileDialog";
            this.ShowInTaskbar = false;
            this.Text = "Profile";
            this.TopMost = true;
            this.FormClosing += new System.Windows.Forms.FormClosingEventHandler(this.ProfileDialog_FormClosing);
            ((System.ComponentModel.ISupportInitialize)(this.chart)).EndInit();
            this.tableLayoutPanel1.ResumeLayout(false);
            this.tableLayoutPanel1.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)(this.numericUpDownProfileWidth)).EndInit();
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.DataVisualization.Charting.Chart chart;
        private System.Windows.Forms.TableLayoutPanel tableLayoutPanel1;
        private System.Windows.Forms.Label labelProfileWidth;
        private System.Windows.Forms.NumericUpDown numericUpDownProfileWidth;
    }
}