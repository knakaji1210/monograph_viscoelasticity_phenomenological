import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

# --- 設定 ---
INITIAL_FREQ = 1.0   
CHANGED_FREQ = 4.0   
CHANGE_TIME = 5.0    
WINDOW_SEC = 3.0     
FPS = 100             
INTERVAL_MS = 1000 / FPS  # 1コマあたりのミリ秒
SAVE_DURATION = 10   
#SAVE_NAME = "./gif/timeLineChart.gif"
SAVE_NAME = "./gif/timeLineChart.mp4"  # mp4で保存する場合

# サイン波の計算関数（絶対時間tから計算）
def calculate_sine_direct(t, freq):
    return np.sin(2 * np.pi * freq * t)

# --- 描画の準備 ---
fig, ax = plt.subplots(figsize=(10, 5))
xdata, ydata = [], []
line, = ax.plot([], [], 'b-', lw=2)

# テキスト表示
time_text = ax.text(0.02, 0.9, '', transform=ax.transAxes, 
                    fontsize=12, fontweight='bold', 
                    bbox=dict(facecolor='white', alpha=0.8))

ax.set_ylim(-1.5, 1.5)
ax.set_xlabel("Time (s)")
ax.grid(True)

DT = 1.0 / FPS

def update(frame):
    current_time = frame * DT
    
    # 周波数の判定
    freq = INITIAL_FREQ if current_time < CHANGE_TIME else CHANGED_FREQ
    
    # 直接計算（位相の積み上げなし）
    y = calculate_sine_direct(current_time, freq)
    
    xdata.append(current_time)
    ydata.append(y)
    
    # メモリ管理その１（表示範囲外をカット）
#    if len(xdata) > int(WINDOW_SEC * 1.5 * FPS):
#        xdata.pop(0)
#        ydata.pop(0)

    # メモリ管理その２
    # 全データを保持せず、現在表示されている範囲（WINDOW_SEC）のデータ＋アルファ
    # だけをリストに残す。これにより、プロット時の計算量が常に一定（低負荷）になる。
    keep_points = int(WINDOW_SEC * FPS * 1.2) 
    if len(xdata) > keep_points:
        del xdata[:-keep_points] # スライスで一気に削除
        del ydata[:-keep_points]

    line.set_data(xdata, ydata)
    time_text.set_text(f"Time: {current_time:.2f}s | Freq: {freq:.1f}Hz")

    # スクロール
    if current_time > WINDOW_SEC:
        ax.set_xlim(current_time - WINDOW_SEC, current_time)
    else:
        ax.set_xlim(0, WINDOW_SEC)
        
    return line, time_text

# --- アニメーションの設定 (intervalを明示的に追加) ---
# frames: 合計フレーム数
# interval: 画面更新間隔（ミリ秒）
ani = FuncAnimation(fig, update, frames=int(SAVE_DURATION * FPS), 
                    blit=True, interval=INTERVAL_MS)

# 保存
print(f"Generating {SAVE_NAME}...")
#ani.save(SAVE_NAME, writer='pillow', fps=FPS)
ani.save(SAVE_NAME, writer='ffmpeg', fps=FPS)  # mp4で保存する場合
print("Finished!")

# 表示確認用
# plt.show()
