# 潮流計算ツール集

MATPOWER互換のNewton-Raphson法による潮流計算ツールです。JSX版とHTML Web版の2つの実装を提供します。

## 🌐 オンラインデモ

**HTML Web版をブラウザで今すぐ試す**: [https://lutelute.github.io/rnw_pf_sample_20251225/powerflow_app.html](https://lutelute.github.io/rnw_pf_sample_20251225/powerflow_app.html)

## 📁 ファイル構成

```
rnw_pf_sample_20251225/
├── power-flow-calculator.jsx    # React JSX版潮流計算ツール
├── powerflow_app.html          # Web HTML版潮流計算ツール
├── case14_init.csv             # IEEE 14母線系統初期値データ
├── case3_init.csv              # 3母線系統初期値データ
├── case9_init.csv              # IEEE 9母線系統初期値データ
└── README.md                   # このファイル
```

## 🚀 機能概要

### React JSX版 (power-flow-calculator.jsx)
- **Newton-Raphson法による高精度潮流計算**
- **MATPOWER互換のデータ形式**
- **リアルタイム系統可視化**
- **IEEE標準テストケース内蔵**

#### 主な機能
- 母線データ編集（PQ/PV/Slack母線対応）
- 送電線データ編集（抵抗・リアクタンス・サセプタンス）
- 発電機データ編集
- Ybus行列表示
- 収束履歴可視化
- 線路潮流・損失計算

### HTML Web版 (powerflow_app.html)
- **ブラウザ単体で動作する完全Web版**
- **時系列シミュレーション機能**
- **配電系統対応**
- **CSV入出力機能**

#### 主な機能
- 豊富なサンプル系統（IEEE 3/9/14/30母線、配電13/33/69ノード）
- インタラクティブ系統図表示
- 時系列負荷変動シミュレーション
- PV発電統合機能
- 電圧プロファイル可視化
- CSV形式でのデータ入出力

## 🛠 使用方法

### JSX版の使用方法
1. React環境を準備
2. `power-flow-calculator.jsx`をプロジェクトにインポート
3. 必要な依存関係をインストール
```bash
npm install react
```

### HTML版の使用方法
1. `powerflow_app.html`をブラウザで開くだけで使用可能
2. サンプル系統を選択して即座に計算開始
3. CSVファイルから独自の系統データをインポート可能

## 📊 サポート系統

### 送電系統
- **IEEE 3母線** - 教育用シンプル系統
- **IEEE 9母線** - WECC 3機9母線系統
- **IEEE 14母線** - 中規模テストケース
- **IEEE 30母線** - 標準テストケース

### 配電系統
- **IEEE 13ノード** - 配電フィーダ
- **IEEE 33母線** - 放射状配電系統
- **IEEE 69母線** - 大規模配電系統

## 💾 データ形式

### CSVファイル形式
プロジェクトに含まれる初期値CSVファイルは、収束困難なケースでの初期値として使用できます。

#### case14_init.csv
```csv
bus,Vm,Va_deg,Pd,Qd,Pg,Qg,type
1,1.060,0.00,0,0,232.4,0,3
2,1.045,-4.98,21.7,12.7,40,0,2
...
```

- `Vm`: 電圧振幅（per unit）
- `Va_deg`: 電圧位相角（度）
- `Pd/Qd`: 有効/無効負荷（MW/MVAr）
- `Pg/Qg`: 有効/無効発電（MW/MVAr）
- `type`: 母線タイプ（1=PQ, 2=PV, 3=Slack）

## ⚡ 計算アルゴリズム

### Newton-Raphson法
- **高速収束**: 通常3-6回で収束
- **高精度**: 1e-6の収束判定
- **ロバスト性**: 複素数演算による数値安定性

### 特徴
1. **複素数ライブラリ**: 独自実装による高速演算
2. **Ybus構築**: π型線路モデル対応
3. **Jacobian行列**: 解析的微分による高精度
4. **ガウス消去法**: 線形方程式求解

## 🎯 対象ユーザー

### 教育
- **電力工学の学習**
- **潮流計算の理解**
- **系統解析の実習**

### 研究開発
- **新アルゴリズム検証**
- **系統計画研究**
- **最適化問題への応用**

### 実務
- **系統運用検討**
- **設備計画**
- **故障時解析**

## 🔧 技術仕様

### JSX版
- **フレームワーク**: React 18+
- **言語**: JavaScript (ES6+)
- **UI**: CSS-in-JS + SVG
- **計算**: 純JavaScript実装

### HTML版
- **技術**: HTML5 + CSS3 + JavaScript
- **可視化**: Canvas API
- **データ処理**: File API
- **UI**: CSS Grid + Flexbox

## ⚠️ 注意事項

### 数値計算
- **単位系**: Per unit基準（Base MVA: 100）
- **収束性**: 初期値に依存する場合あり
- **精度**: 倍精度浮動小数点

### 制限事項
- **タップ比**: 理想変圧器モデル
- **線路**: π型集中定数モデル
- **負荷**: 定インピーダンス

## 📈 拡張可能性

### 将来実装予定
- **最適潮流計算（OPF）**
- **N-1事故時解析**
- **動的シミュレーション**
- **機械学習統合**

## 🤝 貢献

プロジェクトへの貢献を歓迎します：
1. Issues報告
2. Pull Request
3. 機能提案
4. ドキュメント改善

## 📝 ライセンス

このプロジェクトはMITライセンスの下で公開されています。

## 🔗 参考文献

- MATPOWER: https://matpower.org/
- Power System Analysis and Design, Glover et al.
- IEEE Power & Energy Society Test Cases

---

**最終更新**: 2024年12月31日  
**バージョン**: v1.0  
**開発者**: Claude AI Assistant