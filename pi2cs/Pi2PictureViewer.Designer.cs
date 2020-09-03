namespace pi2cs
{
    partial class Pi2PictureViewer
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
            this.panelScrollBarContainer = new System.Windows.Forms.Panel();
            this.scrollBarSlice = new System.Windows.Forms.HScrollBar();
            this.labelSliceNumber = new System.Windows.Forms.Label();
            this.labelPixelValue = new System.Windows.Forms.Label();
            this.panelScrollBarContainer.SuspendLayout();
            this.SuspendLayout();
            // 
            // panelScrollBarContainer
            // 
            this.panelScrollBarContainer.Controls.Add(this.scrollBarSlice);
            this.panelScrollBarContainer.Controls.Add(this.labelSliceNumber);
            this.panelScrollBarContainer.Controls.Add(this.labelPixelValue);
            this.panelScrollBarContainer.Dock = System.Windows.Forms.DockStyle.Bottom;
            this.panelScrollBarContainer.Location = new System.Drawing.Point(0, 410);
            this.panelScrollBarContainer.Margin = new System.Windows.Forms.Padding(4, 4, 4, 4);
            this.panelScrollBarContainer.Name = "panelScrollBarContainer";
            this.panelScrollBarContainer.Size = new System.Drawing.Size(665, 21);
            this.panelScrollBarContainer.TabIndex = 2;
            // 
            // scrollBarSlice
            // 
            this.scrollBarSlice.Dock = System.Windows.Forms.DockStyle.Fill;
            this.scrollBarSlice.Location = new System.Drawing.Point(239, 0);
            this.scrollBarSlice.Name = "scrollBarSlice";
            this.scrollBarSlice.Size = new System.Drawing.Size(426, 21);
            this.scrollBarSlice.TabIndex = 1;
            this.scrollBarSlice.ValueChanged += new System.EventHandler(this.scrollBarSlice_ValueChanged);
            // 
            // labelSliceNumber
            // 
            this.labelSliceNumber.AutoSize = true;
            this.labelSliceNumber.Dock = System.Windows.Forms.DockStyle.Left;
            this.labelSliceNumber.Location = new System.Drawing.Point(149, 0);
            this.labelSliceNumber.Margin = new System.Windows.Forms.Padding(4, 0, 4, 0);
            this.labelSliceNumber.MinimumSize = new System.Drawing.Size(75, 0);
            this.labelSliceNumber.Name = "labelSliceNumber";
            this.labelSliceNumber.Size = new System.Drawing.Size(90, 17);
            this.labelSliceNumber.TabIndex = 2;
            this.labelSliceNumber.Text = "Slice number";
            // 
            // labelPixelValue
            // 
            this.labelPixelValue.AutoSize = true;
            this.labelPixelValue.Dock = System.Windows.Forms.DockStyle.Left;
            this.labelPixelValue.Location = new System.Drawing.Point(0, 0);
            this.labelPixelValue.Margin = new System.Windows.Forms.Padding(4, 0, 4, 0);
            this.labelPixelValue.MinimumSize = new System.Drawing.Size(149, 0);
            this.labelPixelValue.Name = "labelPixelValue";
            this.labelPixelValue.Size = new System.Drawing.Size(149, 17);
            this.labelPixelValue.TabIndex = 0;
            this.labelPixelValue.Text = "(x, y): V";
            // 
            // Pi2PictureViewer
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(8F, 16F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.Controls.Add(this.panelScrollBarContainer);
            this.Margin = new System.Windows.Forms.Padding(4, 4, 4, 4);
            this.Name = "Pi2PictureViewer";
            this.Size = new System.Drawing.Size(665, 431);
            this.KeyPress += new System.Windows.Forms.KeyPressEventHandler(this.Pi2PictureViewer_KeyPress);
            this.panelScrollBarContainer.ResumeLayout(false);
            this.panelScrollBarContainer.PerformLayout();
            this.ResumeLayout(false);

        }

        #endregion
        private System.Windows.Forms.Panel panelScrollBarContainer;
        private System.Windows.Forms.HScrollBar scrollBarSlice;
        private System.Windows.Forms.Label labelSliceNumber;
        private System.Windows.Forms.Label labelPixelValue;
    }
}
