#!/usr/bin/env python3
"""
Simple degradation classifier from pRosettaC ternary complex features.

Uses the best-cluster-per-POI feature table from aggregate_features.py.
With only ~13 samples, fits logistic regression + decision tree + leave-one-out
cross-validation. Reports feature importances and the leave-one-out accuracy.

Usage:
    python degradation_classifier.py
    python degradation_classifier.py --high-affinity-only  # exclude low-affinity POIs
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path

try:
    from sklearn.linear_model import LogisticRegression
    from sklearn.tree import DecisionTreeClassifier, export_text
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import LeaveOneOut
    from sklearn.metrics import accuracy_score, confusion_matrix
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False

PROJECT_ROOT = Path(__file__).resolve().parents[2]
ML_TABLE     = PROJECT_ROOT / 'analysis' / 'degradation_prediction' / 'ml_features_best.csv'
OUT_DIR      = PROJECT_ROOT / 'analysis' / 'degradation_prediction'

# Feature groups for analysis
CONTACT_FEATURES = [
    'best_interface_contacts_5A',
    'best_polar_hbond_pairs',
    'best_interface_area_est_A2',
    'best_poi_iface_atoms',
    'best_vhl_iface_atoms',
]
UBQ_FEATURES = [
    'best_n_lys_e2_optimal',
    'best_e2_lys_dist_1',
    'best_e2_lys_dist_2',
    'md1_best_e2_A',
]
QUALITY_FEATURES = [
    'best_n_poses',
    'n_clusters',
    'best_score_orig_REU',
]
ALL_FEATURES = CONTACT_FEATURES + UBQ_FEATURES + QUALITY_FEATURES

# Core features most likely to discriminate degraders vs non-degraders
CORE_FEATURES = [
    'best_interface_contacts_5A',
    'best_n_lys_e2_optimal',
    'best_polar_hbond_pairs',
    'md1_best_e2_A',
]


def load_table(high_affinity_only=False):
    if not ML_TABLE.exists():
        print(f'Feature table not found: {ML_TABLE}')
        print('Run aggregate_features.py first.')
        return None
    df = pd.read_csv(ML_TABLE)
    if high_affinity_only:
        df = df[df['affinity'] == 'high'].copy()
    return df


def print_feature_summary(df):
    print('\n── Feature summary by label ──')
    for label, grp in df.groupby('degrader'):
        name = 'DEGRADER' if label else 'NON-DEGRADER'
        pois = ', '.join(grp['poi'])
        print(f'\n{name} ({len(grp)}): {pois}')
        for col in CORE_FEATURES:
            if col in grp.columns:
                vals = grp[col].dropna()
                if len(vals):
                    print(f'  {col:40s}: {vals.mean():.1f} ± {vals.std():.1f}  '
                          f'(range {vals.min():.1f}–{vals.max():.1f})')


def loo_evaluate(X, y, poi_names, model, scaler, feature_names):
    """Leave-one-out cross validation."""
    loo  = LeaveOneOut()
    preds, trues, held_out = [], [], []
    for train_idx, test_idx in loo.split(X):
        X_tr, X_te = X[train_idx], X[test_idx]
        y_tr = y[train_idx]
        if scaler:
            sc = StandardScaler()
            X_tr = sc.fit_transform(X_tr)
            X_te = sc.transform(X_te)
        model.fit(X_tr, y_tr)
        preds.append(model.predict(X_te)[0])
        trues.append(y[test_idx[0]])
        held_out.append(poi_names[test_idx[0]])
    acc = accuracy_score(trues, preds)
    cm  = confusion_matrix(trues, preds)
    print(f'  LOO accuracy: {acc:.2f}  ({int(acc*len(y))}/{len(y)} correct)')
    print(f'  Confusion matrix (TN FP / FN TP):\n  {cm}')
    wrong = [(p, int(t), int(pr)) for p, t, pr in zip(held_out, trues, preds) if t != pr]
    if wrong:
        print(f'  Misclassified: {", ".join(f"{p}(true={t},pred={pr})" for p,t,pr in wrong)}')
    return acc, preds, trues


def run_logistic(df, feature_cols, label='degrader'):
    print(f'\n── Logistic Regression (features: {feature_cols}) ──')
    available = [c for c in feature_cols if c in df.columns]
    if not available:
        print('  No matching columns found.')
        return
    sub = df[available + [label, 'poi']].dropna()
    if len(sub) < 4:
        print(f'  Only {len(sub)} samples after dropping NaN — too few.')
        return
    X = sub[available].values
    y = sub[label].values.astype(int)
    pois = sub['poi'].values

    model = LogisticRegression(max_iter=1000, C=0.1)
    loo_evaluate(X, y, pois, model, scaler=True, feature_names=available)

    # Fit on all data for coefficients
    sc = StandardScaler()
    X_sc = sc.fit_transform(X)
    model.fit(X_sc, y)
    print('  Coefficients:')
    for feat, coef in sorted(zip(available, model.coef_[0]), key=lambda x: abs(x[1]), reverse=True):
        print(f'    {feat:40s}: {coef:+.3f}')


def run_decision_tree(df, feature_cols, label='degrader'):
    print(f'\n── Decision Tree (features: {feature_cols}) ──')
    available = [c for c in feature_cols if c in df.columns]
    if not available:
        print('  No matching columns found.')
        return
    sub = df[available + [label, 'poi']].dropna()
    if len(sub) < 4:
        print(f'  Only {len(sub)} samples — too few.')
        return
    X = sub[available].values
    y = sub[label].values.astype(int)
    pois = sub['poi'].values

    model = DecisionTreeClassifier(max_depth=2, random_state=42)
    loo_evaluate(X, y, pois, model, scaler=None, feature_names=available)

    # Fit on all for tree structure
    model.fit(X, y)
    print('\n  Tree structure:')
    print(export_text(model, feature_names=available))
    print('  Feature importances:')
    for feat, imp in sorted(zip(available, model.feature_importances_), key=lambda x: x[1], reverse=True):
        if imp > 0:
            print(f'    {feat:40s}: {imp:.3f}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--high-affinity-only', action='store_true',
                        help='Only use high-affinity POIs (exclude low-affinity degraders)')
    args = parser.parse_args()

    if not HAS_SKLEARN:
        print('scikit-learn not installed. Install with: pip install scikit-learn')
        return

    df = load_table(high_affinity_only=args.high_affinity_only)
    if df is None or df.empty:
        return

    print(f'Loaded {len(df)} POIs ({df["degrader"].sum()} degraders, {(1-df["degrader"]).sum()} non-degraders)')
    print_feature_summary(df)

    # Core model: interface contacts + E2 Lys
    run_logistic(df, CORE_FEATURES)
    run_decision_tree(df, CORE_FEATURES)

    # Contact-only model
    run_logistic(df, CONTACT_FEATURES)

    # UBQ-only model
    run_logistic(df, UBQ_FEATURES)

    # Save summary table
    out = OUT_DIR / 'classification_summary.csv'
    df[['poi', 'degrader', 'affinity'] + [c for c in CORE_FEATURES if c in df.columns]].to_csv(out, index=False)
    print(f'\nSummary table → {out}')


if __name__ == '__main__':
    main()
