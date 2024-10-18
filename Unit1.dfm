object Form1: TForm1
  Left = 306
  Top = 131
  Width = 1057
  Height = 500
  Caption = 'Form1'
  Color = clBtnFace
  Constraints.MinHeight = 500
  Constraints.MinWidth = 750
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  KeyPreview = True
  OldCreateOrder = False
  OnClose = FormClose
  OnCreate = FormCreate
  OnDestroy = FormDestroy
  OnKeyDown = FormKeyDown
  PixelsPerInch = 96
  TextHeight = 13
  object PaintBox1: TPaintBox
    Left = 0
    Top = 0
    Width = 491
    Height = 461
    Cursor = crCross
    Align = alClient
    Color = clBtnFace
    Constraints.MinHeight = 100
    Font.Charset = ANSI_CHARSET
    Font.Color = clBlack
    Font.Height = -13
    Font.Name = 'Consolas'
    Font.Style = []
    ParentColor = False
    ParentFont = False
    OnMouseDown = PaintBox1MouseDown
    OnMouseMove = PaintBox1MouseMove
    OnMouseUp = PaintBox1MouseUp
    OnPaint = PaintBox1Paint
  end
  object Panel1: TPanel
    Left = 491
    Top = 0
    Width = 550
    Height = 461
    Align = alRight
    BevelOuter = bvNone
    TabOrder = 0
    object Memo1: TMemo
      Left = 0
      Top = 89
      Width = 550
      Height = 372
      Align = alClient
      BorderStyle = bsNone
      Color = clBtnFace
      Font.Charset = ANSI_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Consolas'
      Font.Style = []
      Lines.Strings = (
        'Memo1')
      ParentFont = False
      ReadOnly = True
      ScrollBars = ssBoth
      TabOrder = 1
      WordWrap = False
    end
    object Panel2: TPanel
      Left = 0
      Top = 0
      Width = 550
      Height = 89
      Align = alTop
      BevelOuter = bvNone
      TabOrder = 0
      object Label3: TLabel
        Left = 28
        Top = 24
        Width = 53
        Height = 13
        Alignment = taRightJustify
        Caption = 'Frequency '
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = []
        ParentFont = False
      end
      object Label4: TLabel
        Left = 156
        Top = 24
        Width = 25
        Height = 13
        Caption = 'GHz'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = [fsBold]
        ParentFont = False
      end
      object Label6: TLabel
        Left = 36
        Top = 52
        Width = 45
        Height = 13
        Alignment = taRightJustify
        Caption = 'F/D ratio '
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = []
        ParentFont = False
      end
      object Label7: TLabel
        Left = 200
        Top = 52
        Width = 49
        Height = 13
        Alignment = taRightJustify
        Caption = 'Efficiency '
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = []
        ParentFont = False
      end
      object Label8: TLabel
        Left = 373
        Top = 24
        Width = 32
        Height = 13
        Alignment = taRightJustify
        Caption = 'Points '
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = []
        ParentFont = False
      end
      object Label9: TLabel
        Left = 373
        Top = 52
        Width = 32
        Height = 13
        Alignment = taRightJustify
        Caption = 'Petals '
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = []
        ParentFont = False
      end
      object Label11: TLabel
        Left = 324
        Top = 52
        Width = 10
        Height = 13
        Caption = '%'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = [fsBold]
        ParentFont = False
      end
      object Label1: TLabel
        Left = 204
        Top = 24
        Width = 45
        Height = 13
        Alignment = taRightJustify
        Caption = 'Diameter '
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = []
        ParentFont = False
      end
      object Label2: TLabel
        Left = 324
        Top = 24
        Width = 17
        Height = 13
        Caption = 'cm'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = [fsBold]
        ParentFont = False
      end
      object FrequencyEdit: TEdit
        Left = 88
        Top = 20
        Width = 63
        Height = 21
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = [fsBold]
        ParentFont = False
        TabOrder = 0
        Text = '1.69'
        OnKeyDown = EditKeyDown
        OnKeyUp = EditKeyUp
      end
      object FDRatioEdit: TEdit
        Left = 88
        Top = 48
        Width = 63
        Height = 21
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = [fsBold]
        ParentFont = False
        TabOrder = 3
        Text = '0.42'
        OnKeyDown = EditKeyDown
        OnKeyUp = EditKeyUp
      end
      object DiameterEdit: TEdit
        Left = 256
        Top = 20
        Width = 63
        Height = 21
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = [fsBold]
        ParentFont = False
        TabOrder = 1
        Text = '82.5'
        OnKeyDown = EditKeyDown
        OnKeyUp = EditKeyUp
      end
      object PointsCSpinEdit: TCSpinEdit
        Left = 412
        Top = 20
        Width = 61
        Height = 22
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = [fsBold]
        MaxValue = 100
        MinValue = 1
        ParentFont = False
        TabOrder = 2
        Value = 2
        OnChange = CSpinChange
      end
      object PetalsCSpinEdit: TCSpinEdit
        Left = 412
        Top = 48
        Width = 61
        Height = 22
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -11
        Font.Name = 'MS Sans Serif'
        Font.Style = [fsBold]
        MaxValue = 30
        MinValue = 3
        ParentFont = False
        TabOrder = 4
        Value = 12
        OnChange = CSpinChange
      end
    end
    object EfficiencyCSpinEdit: TCSpinEdit
      Left = 256
      Top = 48
      Width = 61
      Height = 22
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -11
      Font.Name = 'MS Sans Serif'
      Font.Style = [fsBold]
      MaxValue = 100
      MinValue = 1
      ParentFont = False
      TabOrder = 2
      Value = 50
      OnChange = CSpinChange
    end
  end
  object OpenDialog1: TOpenDialog
    Filter = 
      'All files|*.*|Supported files *.wav *.raw|*.wav;*.raw|WAV files|' +
      '*.wav|RAW files|*.raw'
    FilterIndex = 2
    Options = [ofReadOnly, ofPathMustExist, ofFileMustExist, ofEnableSizing]
    Left = 20
    Top = 16
  end
  object SaveDialog1: TSaveDialog
    Filter = 
      'All files|*.*|Supported files *.wav *.raw|*.wav;*.raw|WAV files|' +
      '*.wav|RAW files|*.raw'
    FilterIndex = 2
    Options = [ofOverwritePrompt, ofHideReadOnly, ofPathMustExist, ofEnableSizing]
    Left = 20
    Top = 52
  end
end
