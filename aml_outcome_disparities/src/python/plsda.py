"""
    wrap sklearn PLSRegression with ClassifierMixin to get a class that
    can be used for PLS-DA
"""


import numpy as np
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import LabelEncoder


class PLSDAClassifier(BaseEstimator, ClassifierMixin):
    """
    PLS-DA Classifier wrapper around sklearn's PLSRegression.
    Compatible with sklearn's cross_val_score, GridSearchCV, etc.
    """
    
    def __init__(self, n_components=2, scale=True, max_iter=500, tol=1e-06):
        self.n_components = n_components
        self.scale = scale
        self.max_iter = max_iter
        self.tol = tol
    
    def fit(self, X, y):
        # Encode labels
        self.label_encoder_ = LabelEncoder()
        y_encoded = self.label_encoder_.fit_transform(y)
        self.classes_ = self.label_encoder_.classes_
        # For binary: use single column
        # For multiclass: use one-hot encoding
        if len(self.classes_) == 2:
            y_pls = y_encoded
        else:
            # One-hot encode for multiclass
            y_pls = np.zeros((len(y_encoded), len(self.classes_)))  # type: ignore
            for i, val in enumerate(y_encoded):  # type: ignore
                y_pls[i, val] = 1
        # Fit PLS
        self.pls_ = PLSRegression(
            n_components=self.n_components,
            scale=self.scale,
            max_iter=self.max_iter,
            tol=self.tol
        )
        self.pls_.fit(X, y_pls)
        return self
    
    def predict(self, X):
        y_pred_continuous = self.pls_.predict(X)
        if len(self.classes_) == 2:
            # Binary: threshold at 0.5
            y_pred_encoded = (y_pred_continuous.ravel() > 0.5).astype(int)
        else:
            # Multiclass: argmax
            y_pred_encoded = np.argmax(y_pred_continuous, axis=1)
        return self.label_encoder_.inverse_transform(y_pred_encoded)
    
    def predict_proba(self, X):
        """Return pseudo-probabilities (not true probabilities)"""
        y_pred_continuous = self.pls_.predict(X)
        if len(self.classes_) == 2:
            # Convert to 2-column format
            y_pred_continuous = y_pred_continuous.ravel()
            proba = np.column_stack([1 - y_pred_continuous, y_pred_continuous])
        else:
            proba = y_pred_continuous
        # Clip and normalize to [0,1]
        proba = np.clip(proba, 0, 1)
        proba = proba / proba.sum(axis=1, keepdims=True)
        return proba
    
    def transform(self, X):
        """Get PLS scores"""
        return self.pls_.transform(X)
    
    def fit_transform(self, X, y):
        self.fit(X, y)
        return self.transform(X)
    
    @property
    def coef_(self):
        return self.pls_.coef_
    
    @property
    def x_loadings_(self):
        return self.pls_.x_loadings_
    
    @property
    def vip_scores_(self):
        """Calculate VIP (Variable Importance in Projection) scores"""
        t = self.pls_.x_scores_
        w = self.pls_.x_weights_
        q = self.pls_.y_loadings_
        p, h = w.shape
        vips = np.zeros(p)
        s = np.diag(t.T @ t @ q.T @ q).reshape(h, -1)
        total_s = np.sum(s)
        for i in range(p):
            weight = np.array([(w[i, j] / np.linalg.norm(w[:, j]))**2 for j in range(h)])
            vips[i] = np.sqrt(p * (s.T @ weight) / total_s)
        return vips.ravel()