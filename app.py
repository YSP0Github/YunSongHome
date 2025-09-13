from flask import Flask, render_template, request, redirect, url_for, send_from_directory, abort
from flask_sqlalchemy import SQLAlchemy
from datetime import datetime
import os
from werkzeug.utils import secure_filename
from sqlalchemy import or_

# 初始化Flask应用
app = Flask(__name__)
# 配置数据库：当前文件夹下创建 literature.db 文件
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///literature.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False  # 关闭不必要的警告

# 配置文件上传
app.config['UPLOAD_FOLDER'] = os.path.join(app.root_path, 'uploads')  # 上传文件夹路径
app.config['MAX_CONTENT_LENGTH'] = 50 * 1024 * 1024  # 最大50MB
app.config['ALLOWED_EXTENSIONS'] = {'pdf', 'doc', 'docx', 'ppt', 'pptx', 'txt'}  # 允许的文件类型

# 创建上传文件夹（如果不存在）
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
if not os.access(app.config['UPLOAD_FOLDER'], os.W_OK):
    os.chmod(app.config['UPLOAD_FOLDER'], 0o755)
    
db = SQLAlchemy(app)  # 初始化数据库

# 定义文献数据模型（对应数据库表）
class Document(db.Model):
    id = db.Column(db.Integer, primary_key=True)  # 唯一ID
    title = db.Column(db.String(255), nullable=False)  # 文献标题
    authors = db.Column(db.String(255), nullable=False)  # 作者
    journal = db.Column(db.String(255))  # 期刊/会议名
    year = db.Column(db.Integer)  # 发表年份
    volume = db.Column(db.String(50))  # 卷
    issue = db.Column(db.String(50))  # 期
    pages = db.Column(db.String(50))  # 页码
    abstract = db.Column(db.Text)  # 摘要
    category = db.Column(db.String(100))  # 分类（如"计算机科学"）
    tags = db.Column(db.String(255))  # 标签（用逗号分隔）
    view_count = db.Column(db.Integer, default=0)  # 浏览量
    created_at = db.Column(db.DateTime, default=datetime.now)  # 添加时间
    
    # 新增文献文件相关字段
    file_name = db.Column(db.String(255))  # 文件名
    file_path = db.Column(db.String(512))  # 文件存储路径
    file_size = db.Column(db.Integer)  # 文件大小（字节）
    file_type = db.Column(db.String(50))  # 文件类型（如pdf、docx）
    doc_type = db.Column(db.String(50))  # 文献类型（期刊/会议等）

# 检查文件类型是否允许
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in app.config['ALLOWED_EXTENSIONS']

# 生成测试数据（首次运行时执行）
def create_test_data():
    # 先清空现有数据（避免重复添加）
    db.session.query(Document).delete()
    # 添加3条测试文献
    docs = [
        Document(
            title="基于深度学习的图像识别技术在医学诊断中的应用研究",
            authors="张明, 李华, 王芳",
            journal="计算机学报",
            year=2023,
            volume="46",
            issue="8",
            pages="1789-1805",
            abstract="随着深度学习技术的快速发展，图像识别在医学领域的应用越来越广泛...",
            category="计算机科学",
            tags="深度学习,图像识别,医学诊断,CNN",
            view_count=328,
            doc_type="journal"
        ),
        Document(
            title="量子计算在密码学中的安全应用",
            authors="刘强, 赵敏",
            journal="物理学报",
            year=2024,
            volume="73",
            issue="2",
            pages="156-168",
            abstract="量子计算的发展给传统密码学带来了新的挑战与机遇...",
            category="物理学",
            tags="量子计算,密码学,信息安全",
            view_count=156,
            doc_type="journal"
        ),
        Document(
            title="大数据驱动的城市交通流量预测模型",
            authors="陈晓, 杨光, 刘伟",
            journal="人工智能学报",
            year=2022,
            volume="35",
            issue="5",
            pages="890-902",
            abstract="针对城市交通拥堵问题，本文提出一种基于LSTM的交通流量预测模型...",
            category="人工智能",
            tags="大数据,交通预测,LSTM,人工智能",
            view_count=215,
            doc_type="conference"
        )
    ]
    db.session.add_all(docs)
    db.session.commit()

# 首页路由：显示文献列表
@app.route('/')
def index():
    # 从数据库查询所有文献
    documents = Document.query.all()
    # 统计文献总数
    total = len(documents)
    # 把数据传给HTML模板
    return render_template(
        'YSXS.html',
        documents=documents,
        total_count=total
    )

# 文献上传接口
@app.route('/upload', methods=['POST'])
def upload_document():
    # 获取表单数据
    title = request.form.get('doc-title')
    authors = request.form.get('doc-authors')
    journal = request.form.get('doc-journal')
    year = request.form.get('doc-year', type=int)
    volume = request.form.get('doc-volume')
    issue = request.form.get('doc-issue')
    pages = request.form.get('doc-pages')
    category = request.form.get('doc-category')
    doc_type = request.form.get('doc-type')
    keywords = request.form.get('doc-keywords')
    abstract = request.form.get('doc-abstract')
    
    # 验证必填字段
    if not all([title, authors, category, doc_type]):
        return "缺少必填字段", 400
    
    # 处理文件上传
    file = request.files.get('file')
    if not file or file.filename == '':
        return "未选择文件", 400
    
    if file and allowed_file(file.filename):
        # 安全处理文件名（防止路径遍历攻击）
        filename = secure_filename(file.filename)
        # 生成唯一文件名（避免重复）
        unique_filename = f"{datetime.now().strftime('%Y%m%d%H%M%S')}_{filename}"
        file_path = os.path.join(app.config['UPLOAD_FOLDER'], unique_filename)
        
        # 保存文件
        file.save(file_path)
        
        # 创建文献记录
        new_doc = Document(
            title=title,
            authors=authors,
            journal=journal,
            year=year,
            volume=volume,
            issue=issue,
            pages=pages,
            abstract=abstract,
            category=category,
            tags=keywords,
            doc_type=doc_type,
            file_name=filename,
            file_path=file_path,
            file_size=os.path.getsize(file_path),
            file_type=filename.rsplit('.', 1)[1].lower()
        )
        
        db.session.add(new_doc)
        db.session.commit()
        return redirect(url_for('index'))  # 上传成功后返回首页
    
    return "不支持的文件类型", 400

# 检索接口
@app.route('/search')
def search_documents():
    # 获取检索参数
    query = request.args.get('query', '').strip()  # 检索关键词
    category = request.args.get('category', '')    # 分类筛选
    year = request.args.get('year', '')            # 年份筛选
    doc_type = request.args.get('doc_type', '')    # 文献类型筛选
    
    # 基础查询
    search_query = Document.query
    
    # 关键词检索（支持标题、作者、摘要、关键词）
    if query:
        search_query = search_query.filter(
            or_(
                Document.title.like(f'%{query}%'),
                Document.authors.like(f'%{query}%'),
                Document.abstract.like(f'%{query}%'),
                Document.tags.like(f'%{query}%')
            )
        )
    
    # 分类筛选
    if category:
        search_query = search_query.filter(Document.category == category)
    
    # 年份筛选
    if year:
        try:
            search_query = search_query.filter(Document.year == int(year))
        except ValueError:
            pass  # 忽略无效年份格式
    
    # 文献类型筛选
    if doc_type:
        search_query = search_query.filter(Document.doc_type == doc_type)
    
    # 执行查询并统计总数
    documents = search_query.all()
    total = len(documents)
    
    return render_template(
        'YSXS.html',
        documents=documents,
        total_count=total,
        search_query=query  # 回传检索词用于前端显示
    )

# 文献详情页
@app.route('/document/<int:doc_id>')
def document_detail(doc_id):
    # 查询文献详情
    doc = Document.query.get_or_404(doc_id)
    # 增加浏览量
    doc.view_count += 1
    db.session.commit()
    return render_template('document_detail.html', doc=doc)

# 文件下载/在线查看接口
@app.route('/download/<int:doc_id>')
def download_document(doc_id):
    doc = Document.query.get_or_404(doc_id)
    if not os.path.exists(doc.file_path):
        abort(404, description="文件不存在")
    
    # 发送文件（支持在线查看，浏览器会根据文件类型处理）
    return send_from_directory(
        os.path.dirname(doc.file_path),
        os.path.basename(doc.file_path),
        as_attachment=False,  # False表示在线查看，True表示强制下载
        download_name=doc.file_name  # 显示原始文件名
    )

# 启动应用
if __name__ == '__main__':
    # 创建数据库表和测试数据（只在首次运行时执行）
    with app.app_context():
        db.create_all()  # 创建表
        # 检查是否已有数据，没有则添加测试数据
        if Document.query.count() == 0:
            create_test_data()
    # 启动服务器（debug=True表示自动刷新）
    app.run(debug=True)