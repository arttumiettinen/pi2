namespace pi2csWinFormsTest
{
    partial class Form1
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
            this.button1 = new System.Windows.Forms.Button();
            this.labelRetVal = new System.Windows.Forms.Label();
            this.labelErrorMessage = new System.Windows.Forms.Label();
            this.labelCommandList = new System.Windows.Forms.Label();
            this.SuspendLayout();
            // 
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(85, 21);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(75, 23);
            this.button1.TabIndex = 0;
            this.button1.Text = "Run test";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.button1_Click);
            // 
            // labelRetVal
            // 
            this.labelRetVal.AutoSize = true;
            this.labelRetVal.Location = new System.Drawing.Point(43, 75);
            this.labelRetVal.Name = "labelRetVal";
            this.labelRetVal.Size = new System.Drawing.Size(35, 13);
            this.labelRetVal.TabIndex = 1;
            this.labelRetVal.Text = "label1";
            // 
            // labelErrorMessage
            // 
            this.labelErrorMessage.AutoSize = true;
            this.labelErrorMessage.Location = new System.Drawing.Point(371, 75);
            this.labelErrorMessage.Name = "labelErrorMessage";
            this.labelErrorMessage.Size = new System.Drawing.Size(35, 13);
            this.labelErrorMessage.TabIndex = 2;
            this.labelErrorMessage.Text = "label1";
            // 
            // labelCommandList
            // 
            this.labelCommandList.AutoSize = true;
            this.labelCommandList.Location = new System.Drawing.Point(43, 129);
            this.labelCommandList.Name = "labelCommandList";
            this.labelCommandList.Size = new System.Drawing.Size(35, 13);
            this.labelCommandList.TabIndex = 3;
            this.labelCommandList.Text = "label1";
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(628, 816);
            this.Controls.Add(this.labelCommandList);
            this.Controls.Add(this.labelErrorMessage);
            this.Controls.Add(this.labelRetVal);
            this.Controls.Add(this.button1);
            this.Name = "Form1";
            this.Text = "Form1";
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Button button1;
        private System.Windows.Forms.Label labelRetVal;
        private System.Windows.Forms.Label labelErrorMessage;
        private System.Windows.Forms.Label labelCommandList;
    }
}

